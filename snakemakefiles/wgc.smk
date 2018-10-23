import os
import os.path as op
from snakemake.utils import update_config, makedirs
from astropy.table import Table, Column
from snk.utils import get_resource
from snk.utils import make_symlink, check_genome, check_config_default
import yaml

def ndigit(num):
    if num < 1:
        eprint("no digits: %g" % num)
        sys.exit(1)
    digit = 0 
    while num >= 1:
        num /= 10.0
        digit += 1
    return digit

def check_config(c):
    dir_project, dir_cache = c['dir_project'], c['dir_cache']
    c['dirw'] = dir_cache
    dir_raw = "%s/data/raw_output" % dir_project
    if not op.isdir(dir_raw):
        makedirs(dir_raw)
    
    c = check_config_default(c)
    
    dir_cachelink = op.join(dir_project, 'data', 'cache')
    make_symlink(dir_cache, dir_cachelink)
    dir_rawlink = op.join(c['dirw'], c['dird'])
    make_symlink(dir_raw, dir_rawlink)

    c['genomes'] = c['genomes'].split(" ")
    for genome in c['genomes']:
        check_genome(genome, ['blat'], c)
    c['comps'] = [x.split('-') for x in c['comps']]
    return c

configfile: 'config.yaml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    qry = "[a-zA-Z0-9]+",
    tgt = "[a-zA-Z0-9]+",
    idx = "[0-9]+",
    tchrom = "[a-zA-Z0-9]+"

localrules: all, wgc1_prepare, wgc1_break_qry, wgc1_break_tgt

def all_inputs(wildcards):
    inputs = []
    for qry, tgt in config['comps']:
        fo = "%s/%s_%s/23.chain" % (config['wgc']['dir3'], qry, tgt)
        inputs.append(fo)
    return inputs
rule all:
    input:
        all_inputs

include: "rules/wgc.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
