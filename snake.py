#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path as op
import yaml
from astropy.table import Table, Column

def str2bool(v):
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

    for subdir in [c['dirw'], c['variant_dir'], c['temp_dir']]: 
        if not op.isdir(subdir):
            mkdir(subdir)
    #os.chdir(c.dirw)

    for fnv in 'samplelist shared_config'.split():
        fn = c[fnv]
        assert op.isfile(fn), "cannot read %s" % fn

    c['vcf'] = dict()
    c['vbed'] = dict()
    tm = Table(names = ("sid", "gt", "vpre", "opre", "vcf"), dtype = ['O'] * 5)
    t = Table.read(c['samplelist'], format = 'ascii.tab')
    for i in range(len(t)):
        sid, gt = t['sid'][i], t['genotype'][i]
        fv = op.join(c['variant_dir'], "%s.vcf" % gt)
        assert op.isfile(fv), "no vcf found: %s" % fv
        fb = op.join(c['variant_dir'], "%s.bed" % gt)
        assert op.isfile(fb), "no variant-bed found: %s" % fb
        c['vcf'][gt] = fv
        c['vbed'][gt] = fb
    c['t'] = t

    fy = open(c['shared_config'], 'r')
    scfg = yaml.load(fy)
    c.update(scfg)


if __name__ == '__main__':
    print("you lost your way") 
