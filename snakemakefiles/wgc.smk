import os
import os.path as op
from snakemake.utils import update_config, makedirs
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
    c = check_config_default(c)
    c['dirw'] = c['dirc']
    c['comps'] = [x.split('-') for x in c['comps']]
    for genome in set([i for subl in c['comps'] for i in subl]):
        if c['g'][genome]['annotation']:
            check_genome(genome, ['blat','gatk','snpeff'], c)
        else:
            check_genome(genome, ['blat','gatk'], c)
    return c

configfile: 'config.yml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    qry = "[a-zA-Z0-9]+",
    tgt = "[a-zA-Z0-9]+",
    idx = "[0-9]+",
    tchrom = "[a-zA-Z0-9]+"

localrules: all

def all_inputs(wildcards):
    inputs = []
    for qry, tgt in config['comps']:
        dir4 = "%s/%s_%s" % (config['wgc']['dir4'], qry, tgt)
        diro = "%s/%s_%s" % (config['dirr'], qry, tgt)
        inputs.append("%s/10.vnt.bed" % diro)
        inputs.append("%s/15.%s.tsv" % (diro, tgt))
        if config['g'][qry]['annotation']:
            inputs.append("%s/15.%s.tsv" % (diro, qry))
#        inputs.append("%s/10.tsv" % diro)
    return inputs
rule all:
    input:
        all_inputs

include: "rules/wgc.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
