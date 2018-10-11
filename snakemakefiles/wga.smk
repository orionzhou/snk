import os
import os.path as op
from snakemake.utils import update_config, makedirs
from astropy.table import Table, Column
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
    fy = open(c['config_default'], 'r')
    config_default = yaml.load(fy)
    update_config(config_default, c)
    c = config_default
    
    for subdir in [c['dirw'], c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    
    for rsubdir in [c['dirl'], c['dirp']]:
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)
    
    return c

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)

wildcard_constraints:
    qry = "[a-zA-Z0-9]+",
    tgt = "[a-zA-Z0-9]+",
    idx = "[0-9]+",
    tchrom = "[a-zA-Z0-9]+"

def all_inputs(wildcards):
    inputs = []
    for comp in config['comps']:
        qry, tgt = comp['qry'], comp['tgt']
        fo = "%s/%s_%s/23.chain" % (config['wgc']['dirc'], qry, tgt)
        inputs.append(fo)
    return inputs
rule all:
    input:
        all_inputs

include: "rules/wgc.smk"

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))