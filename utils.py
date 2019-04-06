#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import os.path as op
import re
import subprocess as sp
import pandas as pd
import collections
import yaml
from snakemake.utils import update_config, makedirs

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

def get_resource(config, attempt, k1, k2 = ''):
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
    return {'q': q + aq * (attempt - 1),
            'ppn': ppn + appn * (attempt - 1),
            'runtime': runtime + aruntime * (attempt - 1),
            'mem': mem + amem * (attempt - 1),
            'load': load}

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

def check_genome(genome, c):
    dirw = op.join(c['dirg'], genome)
    dbs = list(c['db'].keys())
    dbs = [db for db in dbs if db in c['x'][genome] and c['x'][genome][db]]
    fos = []
    gdic = dict()
    for db in dbs:
        dirx = op.join(dirw, c['db'][db]['xdir'])
        gdic[db] = dict()
        if genome == 'B73' and db == 'hisat2': ### use snp-corrected ref
            dirx = op.join(dirw, '21_dbs/hisat2_snp')
        gdic[db] = dict()
        for k,v in c['db'][db].items():
            if k == 'xdir': continue
            fp = op.join(dirx, v)
            if db == 'snpeff' and k == 'xout': fp = op.join(dirx, genome, v)
            gdic[db][k] = fp

        if db == 'fasta':
            oks = ['ref','chrom_size','chrom_bed','gap']
            if 'region_file' in gdic[db] and op.isfile(gdic[db]['region_file']):
                oks.append('region_file')
                gdic['regions'] = read_region_file(gdic[db]['region_file'])
        elif db == 'annotation':
            oks = ['lgff','lgff_db','rds']
        elif db in ['bowtie2','bwa','bismark','star','hisat2','blastn','blastp','snpeff']:
            oks = ['xout']
        elif db == 'blat':
            oks = ['x.2bit','x.ooc']
        elif db == 'gatk':
            oks = ['xref','xref.dict']
            if 'known_sites' in gdic[db] and op.isfile(gdic[db]['known_sites']):
                oks.append('known_sites')
        else:
            logging.error("unknown db: %s" % db)
            sys.exit(1)
        for ok in oks: fos.append(gdic[db][ok])

    c['g'][genome] = gdic
    for fo in fos:
        assert op.isfile(fo), "%s not found" % fo

def check_config_default(c):
    for fn in [c['config_default']]:
        assert op.isfile(fn), "cannot read %s" % fn

    fy = open(c['config_default'], 'r')
    config_default = yaml.safe_load(fy)
    update_config(config_default, c)
    c = config_default

    dirh0, dirc0, pid = c['dir_project'], c['dir_cache'], c['pid']
    c['dirh'] = op.join(dirh0, pid)
    c['dirc'] = op.join(dirc0, pid)
    c['dird'] = op.join(dirh0, pid, 'data')
    c['dirr'] = op.join(dirc0, pid, 'data')
    dirh, dirc, dird = c['dirh'], c['dirc'], c['dird']
    dir_raw = op.join(c['dird'], 'raw_output')
    for subdir in [dirh,dirc,dird, c['tmpdir'], dir_raw]:
        if not op.isdir(subdir):
            makedirs(subdir)
    make_symlink(dir_raw, c['dirr'])

    dir_cl = op.join(dird, 'cache')
    make_symlink(dirc, dir_cl)

    df = pd.read_excel(c['config_genome'], sheet_name=0, header=0,
        converters={"annotation":bool, "done":bool, "run":bool,
            "fasta":bool, "blat":bool,
            "bwa":bool, "star":bool, "gatk":bool, "hisat2":bool,
            "snpeff":bool, "blastn":bool, "blastp":bool, "bismark":bool})

    xdic, gdic = dict(), dict()
    for i in range(len(df)):
        genome = df['genome'][i]
        xdic[genome] = {x: df[x][i] for x in list(df) if x != 'genome'}
        gdic[genome] = dict()
    c['x'] = xdic
    c['g'] = gdic

    return c

def read_samplelist(samplelist):
    sl = pd.read_csv(samplelist, sep="\t", header=0)
    y1 = dict()
    y1['SampleID'] = sl['SampleID'].tolist()
    y1['t'] = dict()
    y1['gt'] = dict()
    cols = sl.columns.values.tolist()
    for i in range(len(sl)):
        sid = sl['SampleID'][i]
        sdic = {x: sl[x][i] for x in cols}
        y1['t'][sid] = sdic
        if 'Genotype' in sdic and sdic['Genotype']:
            gt = sdic['Genotype']
            if gt not in y1['gt']:
                y1['gt'][gt] = []
            y1['gt'][gt].append(sid)
    y1['Genotypes'] = list(y1['gt'].keys())
    return y1

def check_config_ngs(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']
    for fn in [c['studylist']]:
        assert op.isfile(fn), "cannot read %s" % fn

    df = pd.read_excel(c['studylist'], sheet_name=0, header=0, converters={"meta":bool, "stress":bool, "run": bool})
    y = dict()
    gdic = dict()
    num_run = 0
    for i in range(len(df)):
        #if not df['run'][i]: continue
        if df['run'][i]: num_run += 1
        yid = df['yid'][i]

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], yid, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)
        dir_raw = "%s/08_raw_output/%s" % (c['dird'], yid)
        if not op.isdir(dir_raw): makedirs(dir_raw)
        dir_rawlink = op.join(c['dirw'], yid, 'data')
        make_symlink(dir_raw, dir_rawlink)

        y1 = {x: df[x][i] for x in list(df) if x != 'yid'}
        keys_to_valid = 'source stranded readtype mapper'.split()
        for key_to_valid in keys_to_valid:
            value_to_valid = y1[key_to_valid]
            assert value_to_valid in c['valid'][key_to_valid], "invalid value for key[%s]: %s" % (key_to_valid, value_to_valid)
        fs = "%s/05_read_list/%s.tsv" % (c['dird'], yid)
        assert op.isfile(fs), "samplelist not found: %s" % fs
        y1['samplelist'] = fs
        fsc = "%s/05_read_list/%s.c.tsv" % (c['dird'], yid)
        if not op.isfile(fsc): fsc = fs
        y1['samplelistc'] = fsc
        y1.update(read_samplelist(fs))
        y[yid] = y1

        ref, mapper = y1['reference'], y1['mapper']
        if ref not in gdic:
            gdic[ref] = set()
        gdic[ref].add(mapper)
    c['y'] = y
    #print('working on %s datasets' % len(y.keys()))
    print('working on %s datasets' % num_run)

    for genome in gdic.keys():
       check_genome(genome, c)
    return c

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

if __name__ == '__main__':
    print("you lost your way")
