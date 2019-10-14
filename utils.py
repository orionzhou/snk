#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import os.path as op
import re
import string
from math import ceil
import subprocess as sp
import numpy as np
import pandas as pd
import collections
import yaml
import gspread
from oauth2client.service_account import ServiceAccountCredentials

from snakemake.utils import update_config, makedirs
from jcvi.utils.natsort import natsorted

def ndigit(num):
    if num < 1:
        eprint("no digits: %g" % num)
        sys.exit(1)
    digit = 0
    while num >= 1:
        num /= 10.0
        digit += 1
    return digit

def str2bool(v):
    if isinstance(v, bool):
        return v
    if not isinstance(v, str):
        raise ValueError("invalid literal for boolean: Not a string.")
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    elif v.lower() in ('no', 'false', 'f', '0'):
        return False
    else:
        raise ValueError('invalid literal for boolean: "%s"' % v)

def make_symlink(dst, src):
    if op.islink(src):
        os.system("rm %s" % src)
    elif op.isdir(src):
        os.system("rm -rf %s" % src)
    os.system("ln -s %s %s" % (dst, src))
    if op.islink(op.join(src, op.basename(dst))):
        os.system("rm %s" % op.join(src, op.basename(dst)))

def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n

def get_job_Neo(w, c, rulename):
    N, e, o = (c['job_default'][k] for k in 'N e o'.split())
    ruleid = config[rulename]['id']
    if 'yid' in w:
        if 'sid' in w:
            N = "%s.%s.%s" % (w.yid, ruleid, w.sid)
            e = "%s/%s/%s/%s.e" (w.yid, config['dirj'], ruleid, w.sid)
            o = "%s/%s/%s/%s.o" (w.yid, config['dirj'], ruleid, w.sid)
        elif 'gt' in w:
            if 'rid' in w:
                N = "%s.%s.%s.%s" % (w.yid, ruleid, w.gt, w.rid)
                e = "%s/%s/%s/%s/%s.e" (w.yid, config['dirj'], ruleid, w.gt, w.rid)
                o = "%s/%s/%s/%s/%s.o" (w.yid, config['dirj'], ruleid, w.gt, w.rid)
            else:
                N = "%s.%s.%s" % (w.yid, ruleid, w.gt)
                e = "%s/%s/%s/%s.e" (w.yid, config['dirj'], ruleid, w.gt)
                o = "%s/%s/%s/%s.o" (w.yid, config['dirj'], ruleid, w.gt)
        else:
            N = "%s.%s" % (w.yid, config[rulename]['id'])
            e = "%s/%s/%s.e" (w.yid, config['dirj'], ruleid)
            o = "%s/%s/%s.o" (w.yid, config['dirj'], ruleid)
    elif 'sid' in w:
        N = "%s.%s" % (w.sid, config[rulename]['id'])
        e = "%s/%s/%s.e" (w.sid, config['dirj'], ruleid)
        o = "%s/%s/%s.o" (w.sid, config['dirj'], ruleid)

    return (N, e, o)

def get_task_count(w, c, rulename):
    ntask = None
    if rulename in """fq_dump fastp
        bwa star hisat2
        bam6_markdup bam7a_recal bam7b_bqsr
    """.split():
        assert c.get('y',{}).get(w.yid,{}).get('t',{}).get(w.sid,{}).get('spots') != None, "no spots info for %s" % w.sid
        ntask = c['y'][w.yid]['t'][w.sid]['spots'] / 1000000
    elif rulename in '''
        cv11_hc cv12_merge_vcfs
    '''.split():
        assert c.get('y',{}).get(w.yid,{}).get('gt',{}).get(w.gt) != None, "no Genotype info for %s" % w.gt
        sids = c['y'][w.yid]['gt'][w.gt]
        ntask = sum((c['y'][w.yid]['t'][sid]['spots'] for sid in sids)) / 1000000
    elif rulename in 'bs2_extract'.split():
        assert c.get('y',{}).get(w.yid,{}).get('m',{}).get(w.mid) != None, "no MergeID info for %s" % w.mid
        sids = c['y'][w.yid]['m'][w.mid]
        ntask = sum((c['y'][w.yid]['t'][sid]['spots'] for sid in sids)) / 1000000
    elif rulename in '''
        cv13_rename cv22_genotype_gvcf
    '''.split():
        assert c.get('y',{}).get(w.yid,{}).get('gt',{}).get(w.sid,{}).get('yid') != None, "no yidinfo for %s" % w.sid
        yid0 = config['y'][yid]['t'][sid]['yid']
        gt = config['y'][yid]['t'][sid]['Genotype']
        assert c.get('y',{}).get(yid0,{}).get('gt',{}).get(gt) != None, "no Genotype info for %s" % gt
        sids = c['y'][yid0]['gt'][gt]
        ntask = sum((c['y'][yid0]['t'][sid]['spots'] for sid in sids)) / 1000000
    elif rulename in 'grn2_make'.split():
        assert c.get('t', {}).get(w.nid,{}).get('sample_size') != None, 'no sample size for %s' % w.nid
        ntask = c['t'][w.nid]['sample_size']

    assert ntask != None, 'failed to get task count for %s' % rulename
    return ntask

def get_resource(w, c, rulename):
    ks = 'q ppn runtime mem appn aruntime amem load'.split()
    q, ppn, runtime, mem, appn, aruntime, amem, load = (c['job_default'][k] for k in ks)
    if rulename in c:
        q, ppn, runtime, mem, appn, aruntime, amem, load = (c[rulename][k] for k in ks)

    N, e, o = None, None, None
    if c[rulename]['tasks'] == '':
        return dict(q=q,ppn=ppn,mem=mem,runtime=runtime,appn=appn,amem=amem,aruntime=aruntime,load=load)

    tasks = [int(i) for i in c[rulename]['tasks'].split(',')]
    qs = [i for i in c[rulename]['qs'].split(',')]
    ppns = [int(i) for i in c[rulename]['ppns'].split(',')]
    mems = [int(i) for i in c[rulename]['mems'].split(',')]
    rates = [float(i) for i in c[rulename]['rates'].split(',')]
    n = len(tasks)
    assert n > 1, ">=2 task levels needed"
    if len(qs) == 1: qs = qs * n
    if len(ppns) == 1: ppns = ppns * n
    if len(mems) == 1: mems = mems * n
    if len(rates) == 1: rates = rates * n
    assert len(ppns)==n and len(qs)==n and len(mems)==n and len(rates)==n, "wrong param: qs,ppns,mems,rates"

    ntask = get_task_count(w, c, rulename)
    for i in range(len(tasks)):
        low, high = tasks[i], np.Inf
        if i < len(tasks)-1:
            high = tasks[i+1]

        if ntask >= low and ntask < high:
            q, rate = qs[i], rates[i]
            ppn = max(ppn, ppns[i])
            mem = max(mem, mems[i])
            break

    max_ppn, max_mem, max_runtime = (c['pbs_queue_config'][q][k] for k in 'ppn mem runtime'.split())
    if ppn > max_ppn:
        print('%s: adjusting ppn from %d to %d' % (rulename, ppn, max_ppn))
        ppn = max_ppn
    if mem > max_mem:
        print('%s: adjusting mem from %d to %d' % (rulename, mem, max_mem))
        mem = max_mem

    runtime = max(runtime, ceil(ntask / (ppn * rate)))
    if runtime > max_runtime:
        print('%s: adjusting runtime from %d to %d' % (rulename, runtime, max_runtime))
        runtime = max_runtime

    ppn, mem, runtime = [int(x) for x in (ppn, mem, runtime)]
    return dict(q=q,ppn=ppn,mem=mem,runtime=runtime,appn=appn,amem=amem,aruntime=aruntime,load=load)

def o_get_config_resource(config, k1, k2 = ''):
    assert k1 in config, '%s not in config' % k1
    c = config[k1]
    q, aq = 0, 0
    runtime, aruntime = 8, 2
    mem, amem = 10, 5
    ppn, appn = 1, 0
    load = 1
    if k2 != '':
        assert k2 in config[k1], '%s not in config[%s]' % (k2, k1)
        c2 = config[k1][k2]
        if 'q' in c2: q = c2['q']
        elif 'q' in c: q = c['q']
        if 'aq' in c2: aq = c2['aq']
        elif 'aq' in c: aq = c['aq']
        if 'ppn' in c2: ppn = c2['ppn']
        elif 'ppn' in c: ppn = c['ppn']
        if 'appn' in c2: appn = c2['ppn']
        elif 'appn' in c: appn = c['appn']
        if 'runtime' in c2: runtime = c2['runtime']
        elif 'runtime' in c: runtime = c['runtime']
        if 'aruntime' in c2: aruntime = c2['aruntime']
        elif 'runtime' in c: aruntime = c['aruntime']
        if 'mem' in c2: mem = c2['mem']
        elif 'mem' in c: mem = c['mem']
        if 'amem' in c2: amem = c2['amem']
        elif 'amem' in c: amem = c['amem']
        if 'load' in c2: load = c2['load']
        elif 'load' in c: load = c['load']
    else:
        if 'q' in c: q = c['q']
        if 'aq' in c: aq = c['aq']
        if 'ppn' in c: ppn = c['ppn']
        if 'appn' in c: appn = c['appn']
        if 'runtime' in c: runtime = c['runtime']
        if 'aruntime' in c: aruntime = c['aruntime']
        if 'mem' in c: mem = c['mem']
        if 'amem' in c: amem = c['amem']
        if 'load' in c: load = c['load']
    return dict(q=q, aq=aq, runtime=runtime, aruntime=aruntime,
                mem=mem, amem=amem, ppn=ppn, appn=appn, load=load)

def o_get_resource(w, config, attempt, k1, k2 = ''):
    #r = get_config_resource(config, k1, k2)
    r = config[k1]
    q = r['q'] + r['aq'] * (attempt - 1)
    load = r['load']

    rulename = k1
    if k2 != '': rulename = "_".join((k1, k2))
    r2 = estimate_resource(w, config, rulename)

    ppn = r['ppn']
    if r2.get('ppn') != None:
        ppn = max(ppn, r2['ppn'])

    runtime = r['runtime']
    if r2.get('runtime') != None:
        runtime = max(runtime, r2['runtime'])

    mem = r['mem']
    if r2.get('mem') != None:
        mem = max(mem, r2['mem'])

    ppn = ppn + r['appn'] * (attempt - 1)
    runtime += r['aruntime'] * (attempt - 1)
    mem += r['amem'] * (attempt - 1)

    res = dict(q=q,ppn=ppn,runtime=runtime,mem=mem,load=load)
    for k, v in res.items():
        res[k] = int(v)

    return res

def read_region_file(fr):
    rdic = dict()
    tr = pd.read_csv(fr, sep='\t', header=0)
    chroms = [str(x) for x in range(1,10)]
    for i in range(len(tr)):
        chrom = tr['chrom'][i]
        start = tr['start'][i]
        end = tr['end'][i]
        rid = tr['rid'][i]
        region_str = "%s:%d-%d" % (chrom, start, end)
        rdic[rid] = region_str
#        print("%d regions read for %s" % (len(c[genome]['regions']), genome))
    return rdic

def read_google_sheet(name='jobs', f_cred='/home/springer/zhoux379/.config/google_service_account.json'):
    scope = ['https://spreadsheets.google.com/feeds',
             'https://www.googleapis.com/auth/drive']
    cred = ServiceAccountCredentials.from_json_keyfile_name(f_cred, scope)
    gc = gspread.authorize(cred)

    #book = gc.open_by_key('1SacBnsUW4fzqGYl0k5FVV5AFq2TJvUlBOIWp1NzLh88')
    book= gc.open('coding')
    sheet = book.worksheet(name)
    df = pd.DataFrame(sheet.get_all_records())
    return df

def read_job_config(c):
    ks = 'ppn runtime mem appn aruntime amem load'.split()
    ks2 = 'ppns mems tasks rates'.split()
    cvt = {k: int for k in ks}
    cvt.update({k: str for k in ks2})
    df = pd.read_excel(c['config'], sheet_name='jobs', header=0, converters=cvt)
    df['group'] = df['group'].ffill()
    df = df.fillna(value = c['job_default'])
    jdic = dict()
    for i in range(len(df)):
        rulename= df['rulename'][i]
        jdic[rulename] = {x: df[x][i] for x in list(df) if x != 'rulename'}
        for k in ks:
            jdic[rulename][k] = int(jdic[rulename][k])

    return jdic

def o_read_job_config(c):
    df = read_google_sheet('jobs')
    df = df.replace('', np.nan)
    ks = 'ppn runtime mem appn aruntime amem load'.split()
    ks2 = 'ppns mems tasks rates'.split()
    cvts = {k: int for k in ks}
    cvts.update({k: str for k in ks2})
    df = df.fillna(value = c['job_default'])
    df['group'] = df['group'].ffill()
    df = df.astype(cvts)

    jdic = dict()
    for i in range(len(df)):
        rulename= df['rulename'][i]
        jdic[rulename] = {x: df[x][i] for x in list(df) if x != 'rulename'}
        for k in ks:
            jdic[rulename][k] = int(jdic[rulename][k])

    return jdic

def read_genome_config(c):
    cvts=dict(hybrid=bool, annotation=bool,
        run=bool,
        fasta=bool, blat=bool,
        bwa=bool, star=bool, gatk=bool, hisat2=bool, snpeff=bool,
        lastn=bool, lastp=bool, blastn=bool, blastp=bool,
        bismark=bool,
        tandup=bool, rds=bool)
    df = pd.read_excel(c['config'], sheet_name='genomes', header=0, converters=cvts)
    xdic = dict()
    for i in range(len(df)):
        genome = df['genome'][i]
        xdic[genome] = {x: df[x][i] for x in list(df) if x != 'genome'}

    return xdic

def o_read_genome_config(c):
    df = read_google_sheet('genomes')
    cvts=dict(hybrid=bool, annotation=bool,
        run=bool,
        fasta=bool, blat=bool,
        bwa=bool, star=bool, gatk=bool, hisat2=bool, snpeff=bool,
        lastn=bool, lastp=bool, blastn=bool, blastp=bool,
        bismark=bool,
        tandup=bool, rds=bool)
    df = df.astype(cvts)

    xdic = dict()
    for i in range(len(df)):
        genome = df['genome'][i]
        xdic[genome] = {x: df[x][i] for x in list(df) if x != 'genome'}

    return xdic

def read_study_list(c):
    cvts = dict(interleaved=bool, ase=bool, stress=bool,
                run=bool, runB=bool, runD=bool, runR=bool)
    df = pd.read_excel(c['config'], sheet_name='barn', header=0, converters=cvts)

    defaults = dict(interleaved=False, accession='', format='sra',
                    stranded='no',
                    ase=False, stress=False, ref='Zmays_B73', hisat2='')
    df = df.fillna(defaults)

    for i in range(len(df)):
        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}
        keys_to_valid = 'source format readtype stranded mapper'.split()
        for key_to_valid in keys_to_valid:
            value_to_valid = y1[key_to_valid]
            assert value_to_valid in c['valid'][key_to_valid], "invalid value for key[%s]: %s" % (key_to_valid, value_to_valid)

    return df

def o_read_study_list(c):
    df = read_google_sheet('barn')
    df = df.replace('', np.nan)
    cvts = dict(interleaved=bool, ase=bool, stress=bool,
                done=bool, run=bool, runB=bool, runD=bool, runR=bool)
    defaults = dict(interleaved=False, accession='', format='sra',
                    stranded='no',
                    done=False, run=False, runB=False, runD=False, runR=False,
                    ase=False, stress=False,
                    ref='Zmays_B73', hisat2='')
    df = df.fillna(defaults)
    df = df.astype(cvts)

    for i in range(len(df)):
        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}
        keys_to_valid = 'source format readtype stranded mapper'.split()
        for key_to_valid in keys_to_valid:
            value_to_valid = y1[key_to_valid]
            assert value_to_valid in c['valid'][key_to_valid], "invalid value for key[%s]: %s" % (key_to_valid, value_to_valid)

    return df

def check_genome(genome, c, tag_hisat2=''):
    # if genome not in c['x']:
        # g1, g2 = genome.split('x')
        # if g1 > g2:
            # genome = 'x'.join((g2, g1))
    dirw = op.join(c['dirg'], genome)

    fos = []
    gdic = dict()
    for db, outkeys in c['db']['outkeys'].items():
        if db not in c['x'][genome] or not c['x'][genome][db]:
            continue
        outkeys = outkeys.split()
        xdir = op.join(dirw, c['db'][db]['xdir'])
        if db == 'hisat2' and tag_hisat2.startswith('B73_'):
            xdir = op.join(dirw, '21_dbs/hisat2', tag_hisat2)

        gdic[db] = dict()
        for k,v in c['db'][db].items():
            if k == 'xdir': continue
            fp = op.join(xdir, v)
            if db == 'snpeff' and k == 'xout': fp = op.join(xdir, genome, v)
            gdic[db][k] = fp

        if db == 'fasta':
            for winkey in 'win11 win56 win60 win127'.split():
                fn = gdic[db].get(winkey, '')
                if fn != '':
                    fw = op.join(xdir, fn)
                    if op.isfile(fw):
                        outkeys.append(winkey)
                        gdic[winkey] = read_region_file(fw)
        elif db == 'gatk':
            key = 'known_sites'
            fn = gdic[db].get(key, '')
            if fn != '':
                fw = op.join(xdir, fn)
                if op.isfile(fw):
                    outkeys.append(key)
        elif db == 'annotation':
            key = 'bsseq'
            fn = gdic[db].get(key, '')
            if fn != '':
                fw = op.join(xdir, fn)
                if op.isfile(fw):
                    outkeys.append(key)

        fos += ["%s/%s" % (xdir, c['db'][db][k]) for k in outkeys]

    c['g'][genome] = gdic
    for fo in fos:
        assert op.isfile(fo), "%s not found" % fo

def check_config_default(c):
    for fn in [c['config_default']]:
        assert op.isfile(fn), "cannot read %s" % fn

    cfg_default = yaml.safe_load(open(c['config_default'], 'r'))
    update_config(cfg_default, c)
    c = cfg_default

    for fn in [c['config_job_default']]:
        assert op.isfile(fn), "cannot read %s" % fn

    cfg_job = read_job_config(c)
    update_config(cfg_job, c)
    c = cfg_job

    dirh0, dirc0 = c['dir_project'], c['dir_cache']
    pid, wid, oid = c['pid'], c['wid'], c['oid']
    c['dirh'] = op.join(dirh0, pid, wid)
    c['dirc'] = op.join(dirc0, pid, wid)
    c['dirr'] = op.join(dirh0, pid, wid, oid)
    dirh, dirc, dirr = c['dirh'], c['dirc'], c['dirr']
    dirr_l = op.join(dirc, oid)
    for subdir in [dirh,dirc,dirr, c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    make_symlink(dirr, dirr_l)

    dirh_l = op.join(dirc, 'primary')
    dirc_l = op.join(dirh, 'cache')
    make_symlink(dirc, dirc_l)
    make_symlink(dirh, dirh_l)

    xdic = read_genome_config(c)
    gdic = {g: dict() for g in xdic.keys()}
    c['x'] = xdic
    c['g'] = gdic

    return c

def read_samplelist(diri, yid, part_size=100000000, cap_gt=False):
    y1 = dict()
    fs = "%s/%s.tsv" % (diri, yid)
    #fsc = "%s/%s.c.tsv" % (diri, yid)
    assert op.isfile(fs), "samplelist not found: %s" % fs
    y1['samplelist'] = fs
    #if not op.isfile(fsc): fsc = fs
    #y1['samplelistc'] = fsc

    sl = pd.read_csv(fs, sep="\t", header=0)
    y1['SampleID'] = sl['SampleID'].tolist()
    y1['t'], y1['gt'], y1['m'] = dict(), dict(), dict()
    cols = sl.columns.values.tolist()
    for i in range(len(sl)):
        sid = sl['SampleID'][i]
        sdic = {x: sl[x][i] for x in cols}
        npart = ceil(sdic['spots'] / part_size)
        balanced_part_size = ceil(sdic['spots'] / npart)
        if npart > 26:
            print("too many pieces (>26) to split: %s|%s" % (yid, sid))
            sys.exit(1)
        sdic['npart'] = npart
        sdic['part_size'] = balanced_part_size
        sdic['parts'] = list(string.ascii_lowercase)[:npart]
        sdic['p'] = {sdic['parts'][i]: i+1 for i in range(sdic['npart'])}
        if cap_gt:
            sdic['Genotype'] = sdic['Genotype'].upper().replace('X', 'x')
        y1['t'][sid] = sdic

        if 'Genotype' in sdic and sdic['Genotype']:
            gt = sdic['Genotype']
            if gt not in y1['gt']:
                y1['gt'][gt] = []
            y1['gt'][gt].append(sid)

        if 'MergeID' in sdic and sdic['MergeID']:
            mid = sdic['MergeID']
            if mid not in y1['m']:
                y1['m'][mid] = []
            y1['m'][mid].append(sid)

    y1['Genotypes'] = list(y1['gt'].keys())
    y1['MergeID'] = list(y1['m'].keys())
    return y1

def check_config_barn(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    df = read_study_list(c)
    y, num_run = dict(), 0
    for i in range(len(df)):
        if not df['run'][i]: continue
        if df['run'][i]: num_run += 1
        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}
        yid = df['yid'][i]

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        idir = c['barn']['idir_sra']
        if y1['source'] == 'local':
            idir = c['barn']['idir_local']
        fs = "%s/%s/%s.tsv" % (c['dirh'], idir, yid)
        assert op.isfile(fs), "samplelist not found: %s" % fs
        y1['samplelist'] = fs

        sl = pd.read_csv(fs, sep="\t", header=0)
        y1['SampleID'] = sl['SampleID'].tolist()
        y1['t'] = dict()
        cols = sl.columns.values.tolist()
        for i in range(len(sl)):
            sid = sl['SampleID'][i]
            sdic = {x: sl[x][i] for x in cols}
            y1['t'][sid] = sdic

        y[yid] = y1

    c['y'] = y
    print('working on %s datasets' % num_run)

    return c

def check_config_rnaseq(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    df = read_study_list(c)
    y, refs = dict(), set()
    num_run = 0
    for i in range(len(df)):
        if not df['runR'][i]: continue
        if df['runR'][i]: num_run += 1
        yid = df['yid'][i]
        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        diri = op.join(c['barn']['home'], c['barn']['odir'])
        y1.update(read_samplelist(diri = diri, yid = yid, cap_gt = True))
        y[yid] = y1

        ref, tag_hisat2 = y1['ref'], y1['hisat2']
        if ref not in refs:
            check_genome(ref, c, tag_hisat2)
        refs.add(y1['ref'])

    c['y'] = y
    print('working on %s datasets' % num_run)

    return c

def create_sample_id_mapping(gts, yid, fo):
    fho = open(fo, 'w')
    print("writing sample ID mapping for %d genotypes" % len(gts))
    for gt in gts:
        ngt = "%s|%s" % (yid, gt)
        print("%s %s" % (gt, ngt), file = fho)
    fho.close()

def check_config_reseq(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    df = read_study_list(c)
    y, refs = dict(), set()
    num_run = 0
    for i in range(len(df)):
        if not df['runD'][i]: continue
        if df['runD'][i]: num_run += 1
        yid = df['yid'][i]
        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        diri = op.join(c['barn']['home'], c['barn']['odir'])
        y1.update(read_samplelist(diri = diri, yid = yid))
        y[yid] = y1

        refs.add(y1['ref'])

    c['y'] = y
    print('working on %s datasets' % num_run)

    for ref in refs:
       check_genome(ref, c)
    return c

def check_config_reseq2(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']
    y, refs = dict(), set()
    num_run = 0
    for yid, ydic in c['y'].items():
        if not ydic['run']: continue
        num_run += 1

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        ydic['t'] = dict()
        fs = "%s/%s/%s.tsv" % (c['dirh'], c['callvnt2']['idir'], yid)
        sl = pd.read_csv(fs, sep="\t", header=0)
        cols = sl.columns.values.tolist()
        for i in range(len(sl)):
            sid = sl['sid'][i]
            sdic = {x: sl[x][i] for x in cols if x != 'sid'}
            ydic['t'][sid] = sdic

        refs.add(ydic['ref'])
        y[yid] = ydic

    for ref in refs:
       check_genome(ref, c)

    c['y'] = y
    print('working on %s datasets' % num_run)

    return c

def check_config_reseq3(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']
    y, refs = dict(), set()
    num_run = 0
    for yid, ydic in c['y'].items():
        if not ydic['run']: continue
        num_run += 1

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        fs = "%s/31_vnt_list/%s.txt" % (c['dirh'], yid)
        assert op.isfile(fs), "samplelist not found: %s" % fs
        sl = pd.read_csv(fs, sep="\t", names=['sid'], header=None)
        ydic['samples'] = sl['sid'].tolist()

        fs = "%s/35_vnt_ase/%s.tsv" % (c['dirh'], yid)
        sl = pd.read_csv(fs, sep="\t", header=0)
        ydic['ase_genotypes'] = sl['Genotype'].tolist()
        cols = sl.columns.values.tolist()
        ydic['ase'] = dict()
        for i in range(len(sl)):
            gt = sl['Genotype'][i]
            gdic = {x: sl[x][i] for x in cols if x != 'Genotype'}
            ydic['ase'][gt] = gdic

        vid = ydic['vid']
        fv = "/home/springer/zhoux379/projects/reseq/data/cache/vcf/%s.vcf.gz" % vid
        assert op.isfile(fv), "vcf not found: %s" % fv
        ydic['vcf'] = fv

        refs.add(ydic['ref'])
        y[yid] = ydic

    for ref in refs:
       check_genome(ref, c)

    c['y'] = y
    print('working on %s datasets' % num_run)

    return c

def check_config_bsseq(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    df = read_study_list(c)
    y, refs = dict(), set()
    num_run = 0
    meta = False
    for i in range(len(df)):
        if not df['runB'][i]: continue
        if df['runB'][i]: num_run += 1
        yid = df['yid'][i]
        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        diri = op.join(c['barn']['home'], c['barn']['odir'])
        y1.update(read_samplelist(diri=diri, yid=yid, part_size=c['trimming']['part_size']))
        y[yid] = y1

        refs.add(y1['ref'])

    c['y'] = y
    print('working on %s datasets' % num_run)

    for ref in refs:
       check_genome(ref, c)
    return c

def check_config_phylo(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']
    for fn in [c['studylist']]:
        assert op.isfile(fn), "cannot read %s" % fn

    df = pd.read_excel(c['studylist'], sheet_name=0, header=0,
                       converters={"run":bool})

    y, refs = dict(), set()
    num_run = 0
    for yid, ydic in c['y'].items():
        if not ydic['run']: continue
        num_run += 1

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

        fs = "%s/05_sample_list/%s.txt" % (c['dirh'], yid)
        assert op.isfile(fs), "samplelist not found: %s" % fs
        sl = pd.read_csv(fs, sep="\t", names=['sid'], header=None)
        ydic['samples'] = sl['sid'].tolist()

        vid = ydic['vid']
        fv = "/home/springer/zhoux379/projects/reseq/data/cache/vcf/%s.vcf.gz" % vid
        assert op.isfile(fv), "vcf not found: %s" % fv
        ydic['vcf'] = fv

        refs.add(ydic['ref'])
        y[yid] = ydic

    for ref in refs:
       check_genome(ref, c)

    c['y'] = y
    print('working on %s datasets' % num_run)

    return c

def check_config_popgen(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    fs = c['popgen']['sample_list']
    sids = [line.rstrip('\n') for line in open(fs,'r')]
    c['popgen']['sids'] = sids

    fg = c['popgen']['gene_id']
    gids = [line.rstrip('\n') for line in open(fg,'r')]
    n_batch = ceil(len(gids) / c['popgen']['batch_size'])
    c['popgen']['batches'] = list(range(1,n_batch+1))

    print("working on %s samples" % len(sids))
    print("processing %s genes in %d batches" % (len(gids), n_batch))
    return c

def check_config_grn(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    f_cfg = op.join(c['dirh'], c['grn']['cfg'])
    df = pd.read_excel(f_cfg, sheet_name=0, header=0)
    c['nid'] = []
    c['t'] = dict()
    for i in range(len(df)):
        nid = df['nid'][i]
        if not nid.startswith('np_gibeerish'):
            c['nid'].append(nid)
            c['t'][nid] = {x: df[x][i] for x in list(df)}
    return c

def check_config_wgc(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']
    c['comps'] = [x.split('-') for x in c['comps']]
    genomes = [i for subl in c['comps'] for i in subl]
    for genome in set(genomes + c['ortho_genomes']):
        check_genome(genome, c)
    return c

if __name__ == '__main__':
    print("you lost your way")
