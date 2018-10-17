#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import os.path as op
import re
import subprocess as sp
import collections
import yaml
from snakemake.utils import update_config, makedirs
from astropy.table import Table, Column

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
    if op.isdir(src):
        os.system("rm -rf %s" % src)
    elif op.islink(src):
        os.system("rm %s" % src)
    os.system("ln -sf %s %s" % (dst, src))

def check_genome(genome, dbs, c):
    dirw = op.join(c['dirg'], c[genome]['gdir'])
    ref, size, regions, gtf = [op.join(dirw, x) for x in [
        '10_genome.fna',
        '15_intervals/01.chrom.sizes',
        '15_intervals/20.gap.sep.60win.tsv',
        '50_annotation/10.gtf']]
    for fi in [ref, size, gtf]:
        assert op.isfile(fi), "%s not found" % fi
    c[genome]['ref'] = ref
    c[genome]['size'] = size 
    c[genome]['gtf'] = gtf
    if op.isfile(regions):
        fr = regions
        c[genome]['regions'] = dict() 
        tr = Table.read(fr, format = 'ascii.tab')
        chroms = [str(x) for x in range(1,10)]
        for i in range(len(tr)):
            chrom = tr['chrom'][i]
            start = tr['start'][i]
            end = tr['end'][i]
            rid = tr['rid'][i]
            region_str = "%s:%d-%d" % (chrom, start, end)
            c[genome]['regions'][rid] = region_str
        print("%d regions read for %s" % (len(c[genome]['regions']), genome))

    if isinstance(dbs, str): dbs = [dbs]
    if 'bwa' in dbs and 'gatk' not in dbs:
        dbs.append('gatk')
    for db in dbs:
        dirx = op.join(dirw, '21_dbs', c[db]['xdir'])
        fos = []
        if db in ['star', 'bwa', 'hisat2']:
            xpre, xout = c[db]['xpre'], c[db]['xout']
            fp, fo = op.join(dirx, xpre), op.join(dirx, xout)
            if genome == 'B73' and db == 'hisat2': ### use snp-corrected ref
                dirx1 = op.join(dirw, '21_dbs', 'hisat2_snp')
                fp, fo = op.join(dirx1, xpre), op.join(dirx1, xout)
            c[genome][db] = fp
            fos = [fo]
        elif db == 'blat':
            ks = ['x.2bit', 'x.ooc']
            c[genome][db] = {x: op.join(dirx, c[db][x]) for x in ks}
            fos = list(c[genome][db].values())
        elif db == 'gatk':
            ks = ['xref', 'xref.fai', 'xref.dict']
            c[genome][db] = {x: op.join(dirx, c[db][x]) for x in ks}
            fos = list(c[genome][db].values())
        else:
            print("unknown db: %s" % db)
            sys.exit(1)
        for fo in fos:
            assert op.isfile(fo), "%s not found" % fo

def check_config(c):
    for fn in [c['studylist'], c['config_default']]:
        assert op.isfile(fn), "cannot read %s" % fn

    fy = open(c['config_default'], 'r')
    config_default = yaml.load(fy)
    update_config(config_default, c)
    c = config_default
    
    study = c['study']
    t = Table.read(c['studylist'], format = 'ascii.tab')
    dic_study = { t['sid'][i]: {x: t[x][i] for x in t.colnames} for i in range(len(t)) }
    assert study in dic_study, "study not in config file: %s" % study
    c['source'] = dic_study[study]['source']
    c['readtype'] = dic_study[study]['readtype']
    c['stranded'] = dic_study[study]['stranded']
    c['reference'] = dic_study[study]['reference']
    c['mapper'] = dic_study[study]['mapper']
    assert c['source'] in ['sra', 'local', 'local_interleaved'], "unknown source: %s" % c['source']
    assert c['stranded'] in ['yes', 'no', 'reverse'], "unknown strand: %s" % c['stranded']
    assert c['readtype'] in ['illumina', 'solid', '3rnaseq'], "unknown readtype: %s" % c['readtype']
    assert c['mapper'] in ['star', 'hisat2', 'bwa'], "unknown mapper: %s" % c['mapper']
    c['genome'] = c['reference']
    dbs = [c['mapper']]
    check_genome(c['reference'], dbs, c)
   
    dir_project, dir_cache = c['dir_project'], c['dir_cache']
    dirw = op.join(dir_cache, study)
    samplelist = "%s/data/05_read_list/%s.tsv" % (dir_project, study)
    assert op.isfile(samplelist), "samplelist not found: %s" % samplelist
    dir_raw = "%s/data/08_raw_output/%s" % (dir_project, study)
    c['dirw'], c['samplelist'] = dirw, samplelist
   
    for subdir in [c['dirw'], dir_raw, c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    for rsubdir in [c['dirl'], c['dirp']]: 
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)
    
    dir_cachelink = op.join(dir_project, 'data', 'cache')
    make_symlink(dir_cache, dir_cachelink)
    dir_rawlink = op.join(dirw, c['dird'])
    make_symlink(dir_raw, dir_rawlink)

    t = Table.read(c['samplelist'], format = 'ascii.tab')
    #tm = Table(names = ("sid", "gt", "vpre", "opre", "vcf"), dtype = ['O'] * 5)
    c['SampleID'] = t['SampleID']
    c['t'] = dict()
    cols = t.colnames
    for i in range(len(t)):
        sid = t['SampleID'][i]
        sdic = {x: t[x][i] for x in cols}
        if 'paired' in sdic:
            sdic['paired'] = str2bool(sdic['paired'])
        c['t'][sid] = sdic

    return c

def check_config_ase(c):
    for subdir in [c['ase']['vdir']]: 
        if not op.isdir(subdir):
            mkdir(subdir)
    t = c['t']
    c['vcf'] = dict()
    c['vbed'] = dict()
    for sid in c['SampleID']:
        gt = t[sid]['Genotype']
        fv = op.join(c['ase']['vdir'], "%s.vcf" % gt)
        fb = op.join(c['ase']['vdir'], "%s.bed" % gt)
        if not op.isfile(fb):
            fv = op.join(c['ase']['vdir2'], "%s.vcf" % gt)
            fb = op.join(c['ase']['vdir2'], "%s.bed" % gt)
        assert op.isfile(fv), "no vcf found: %s" % fv
        assert op.isfile(fb), "no variant-bed found: %s" % fb
        c['vcf'][sid] = fv
        c['vbed'][sid] = fb

# From https://github.com/giampaolo/psutil/blob/master/scripts/meminfo.py
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
