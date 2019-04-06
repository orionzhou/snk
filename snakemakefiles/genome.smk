import os
import os.path as op
import pandas as pd
from snakemake.utils import update_config, makedirs
from snk.utils import get_resource, check_config_default
import yaml

def check_config(c):
    c = check_config_default(c)
    c['dirw'] = c['dird']

    df = pd.read_excel(c['config_genome'], sheet_name=0, header=0,
        converters={"annotation":bool, "done":bool, "run":bool,
            "fasta":bool, "blat":bool,
            "bwa":bool, "star":bool, "gatk":bool, "hisat2":bool,
            "snpeff":bool, "blastn":bool, "blastp":bool, "bismark":bool})

    num_run = 0
    gdic = dict()
    for i in range(len(df)):
        if df['run'][i]: num_run += 1
        genome = df['genome'][i]

        for rsubdir in [c['dirl'], c['dirj']]:
            subdir = op.join(c['dirw'], genome, rsubdir)
            if not op.isdir(subdir):
                makedirs(subdir)
        gdic[genome] = {x: df[x][i] for x in list(df) if x != 'genome'}
    c['x'] = gdic
    print('working on %s genomes' % num_run)

    return c

configfile: 'config.yml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    genome = "[a-zA-Z0-9]+",
    chrom = "[a-zA-Z0-9]+"

def all_inputs(wildcards):
    inputs = []
    for genome in config['x']:
        if not config['x'][genome]['run']: continue
        dbs = ['fasta','annotation','blat','bwa','star','gatk','hisat2','snpeff',
            'blastn','blastp','bismark']
        dbs = [db for db in dbs if config['x'][genome][db]]
        for db in dbs:
            odir = "%s/%s" % (genome, config['db'][db]['xdir'])
            if db == 'fasta':
                ks = ['ref','chrom_size','chrom_bed','gap','fchain','bchain']
            elif db == 'annotation':
                ks = ['gff','gff_db','gtf','tsv','bed','des','fna','faa',
                    'lgff','lgff_db','lgtf','ltsv','lbed','ldes','lfna','lfaa',
                    'tandup','rds']
            elif db in ['bowtie2','bwa','bismark','star','hisat2','blastn','blastp']:
                ks = ['xout']
            elif db == 'blat':
                ks = ['x.2bit','x.ooc']
            elif db == 'gatk':
                ks = ['xref','xref.dict']
            elif db == 'snpeff':
                ks = ['xcfg']
                inputs.append("%s/%s/%s/%s" % (genome, config['db'][db]['xdir'], genome, config['db'][db]['xout']))
            inputs += ["%s/%s/%s" % (genome, config['db'][db]['xdir'], config['db'][db][k]) for k in ks]
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
