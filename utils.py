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
    c['paired'] = str2bool(c['paired'])
    assert c['stranded'] in ['yes', 'no', 'reverse'], "unknown strand option: %s" % c['stranded']

    fy = open(c['config_default'], 'r')
    config_default = yaml.load(fy)
    update_config(config_default, c)
    c = config_default
    
    for subdir in [c['dirw'], c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    
    for rsubdir in [c['dirl'], c['dirp'], c['dirq']]: 
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)

    for fn in [c['samplelist'], c['config_default']]:
        assert op.isfile(fn), "cannot read %s" % fn

    tm = Table(names = ("sid", "gt", "vpre", "opre", "vcf"), dtype = ['O'] * 5)
    t = Table.read(c['samplelist'], format = 'ascii.tab')
    c['t'] = t

    return c

def check_config_ase(c):
    for subdir in [c['ase']['variant_dir']]: 
        if not op.isdir(subdir):
            mkdir(subdir)
    t = c['t']
    c['vcf'] = dict()
    c['vbed'] = dict()
    for i in range(len(t)):
        sid, gt = t['sid'][i], t['genotype'][i]
        fv = op.join(c['ase']['variant_dir'], "%s.vcf" % gt)
        fb = op.join(c['ase']['variant_dir'], "%s.bed" % gt)
        if not op.isfile(fb):
            fv = op.join(c['ase']['variant_dir2'], "%s.vcf" % gt)
            fb = op.join(c['ase']['variant_dir2'], "%s.bed" % gt)
        #assert op.isfile(fv), "no vcf found: %s" % fv
        assert op.isfile(fb), "no variant-bed found: %s" % fb
        #c['vcf'][gt] = fv
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
