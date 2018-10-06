import os
import os.path as op
from snakemake.utils import update_config, makedirs
from astropy.table import Table, Column
import yaml

def check_config(c):
    fy = open(c['config_default'], 'r')
    config_default = yaml.load(fy)
    update_config(config_default, c)
    c = config_default
    
    for subdir in [c['dirw'], c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    
    #for rsubdir in [c['dirl'], c['dirp']]:
    for rsubdir in [c['dirp']]:
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)
    
    return c

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)

wildcard_constraints:
    genome = "[a-zA-Z0-9]+",
    chrom = "[a-zA-Z0-9]+"

def all_inputs(wildcards):
    inputs = []
    for genome in config['genome']:
        dbs = set(list(config['genome'][genome]['dbs'].split()))
        for db in dbs:
            odir = genome if db == 'fasta' else "%s/21_dbs/%s" % (genome, config[db]['odir'])
            if db == 'fasta':
                odir = "%s" % genome
                inputs.append("%s/10_genome.fna" % odir)
                inputs.append("%s/10_genome.fna.fai" % odir)
                inputs.append("%s/15_intervals/01.chrom.bed" % odir)
                inputs.append("%s/15_intervals/01.chrom.sizes" % odir)
                inputs.append("%s/15_intervals/11.gap.bed" % odir)
            if db == 'blat':
                inputs.append("%s/db.2bit" % odir)
            if db == 'bowtie2':
                inputs.append("%sdb.rev.1.bt2" % odir)
            if db == 'bwa':
                inputs.append("%s/db.bwt" % odir)
            if db == 'star':
                inputs.append("%s/SA" % odir)
            if db == 'gatk':
                inputs.append("%s/db.fasta" % odir)
                inputs.append("%s/db.dict" % odir)
    return inputs
rule all:
    input:
        all_inputs

include: "rules/genome.smk"

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
