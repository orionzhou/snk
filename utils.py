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

def check_config(c):
    assert c['stranded'] in ['yes', 'no', 'reverse'], "unknown strand option: %s" % c['stranded']

    fy = open(c['config_default'], 'r')
    config_default = yaml.load(fy)
    update_config(config_default, c)
    c = config_default
    
    for subdir in [c['dirw'], c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    
    for rsubdir in [c['dirl'], c['dirp'], c['dird']]: 
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)

    for fn in [c['samplelist'], c['config_default']]:
        assert op.isfile(fn), "cannot read %s" % fn

    tm = Table(names = ("sid", "gt", "vpre", "opre", "vcf"), dtype = ['O'] * 5)
    t = Table.read(c['samplelist'], format = 'ascii.tab')
    c['SampleID'] = t['SampleID']
    c['t'] = dict()
    cols = t.colnames
    
    for i in range(len(t)):
        sid = t['SampleID'][i]
        sdic = {x: t[x][i] for x in cols}
        if 'paired' in sdic:
            sdic['paired'] = str2bool(sdic['paired'])
        c['t'][sid] = sdic

    if 'regions' in c['db'] and op.isfile(c['db']['regions']):
        fr = c['db']['regions']
        c['regions'] = dict()
        tr = Table.read(fr, format = 'ascii.tab')
        chroms = [str(x) for x in range(1,10)]
        for i in range(len(tr)):
            chrom = tr['chrom'][i]
            start = tr['start'][i]
            end = tr['end'][i]
            rid = tr['rid'][i]
            region_str = "%d:%d-%d" % (chrom, start, end)
            c['regions'][rid] = region_str
        print("%d regions read" % len(c['regions'].keys()))

    return c

def check_config_ase(c):
    for subdir in [c['ase']['variant_dir']]: 
        if not op.isdir(subdir):
            mkdir(subdir)
    t = c['t']
    c['vcf'] = dict()
    c['vbed'] = dict()
    for sid in c['SampleID']:
        gt = t[sid]['Genotype']
        fv = op.join(c['ase']['variant_dir'], "%s.vcf" % gt)
        fb = op.join(c['ase']['variant_dir'], "%s.bed" % gt)
        if not op.isfile(fb):
            fv = op.join(c['ase']['variant_dir2'], "%s.vcf" % gt)
            fb = op.join(c['ase']['variant_dir2'], "%s.bed" % gt)
        assert op.isfile(fv), "no vcf found: %s" % fv
        assert op.isfile(fb), "no variant-bed found: %s" % fb
        c['vcf'][gt] = fv
        c['vbed'][gt] = fb

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
