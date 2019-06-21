import os
import os.path as op
import pandas as pd
from snakemake.utils import update_config, makedirs
from snk.utils import get_resource, check_config_default, read_genome_config, make_symlink
import yaml

def check_config(c):
    c = check_config_default(c)
    c['dirw'] = c['dirh']
    xdic = read_genome_config(c)

    num_run = 0
    for genome in xdic.keys():
        if xdic[genome]['run']:
            num_run += 1
        else:
            continue

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], genome, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)

    c['x'] = xdic
    print('working on %s genomes' % num_run)
    return c

configfile: 'config.yml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    genome = "[a-zA-Z0-9_]+",
    chrom = "[a-zA-Z0-9]+"

def all_inputs(wildcards):
    inputs = []
    for genome in config['x']:
        if not config['x'][genome]['run']:
            continue
        for db, outkeys in config['db']['outkeys'].items():
            if db not in config['x'][genome] or not config['x'][genome][db]:
                continue
            odir = "%s/%s" % (genome, config['db'][db]['xdir'])
            outkeys = outkeys.split()
            inputs += ["%s/%s" % (odir, config['db'][db][k]) for k in outkeys]
            if db == 'snpeff':
                inputs.append("%s/%s/%s" % (odir, genome, config['db'][db]['xout']))

    return inputs
localrules: all
rule all:
    input:
        all_inputs

include: "rules/genome.seq.smk"
include: "rules/genome.gff.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
