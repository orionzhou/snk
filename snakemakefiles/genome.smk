import os
import os.path as op
from snakemake.utils import update_config, makedirs
from astropy.table import Table, Column
import yaml
from snk.utils import get_resource

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
    for genome in config['genomes']:
        dbs = set(list(config[genome]['dbs'].split()))
        for db in dbs:
            odir = genome if db == 'fasta' else "%s/21_dbs/%s" % (genome, config[db]['xdir'])
            if db == 'fasta':
                odir = "%s" % genome
                inputs.append("%s/10_genome.fna" % odir)
                inputs.append("%s/10_genome.fna.fai" % odir)
                inputs.append("%s/15_intervals/01.chrom.bed" % odir)
                inputs.append("%s/15_intervals/01.chrom.sizes" % odir)
                inputs.append("%s/15_intervals/11.gap.bed" % odir)
            if db in ['bowtie2','bwa','bismark','star','hisat2','blastn','blastp']:
                inputs.append("%s/%s" % (odir, config[db]['xout']))
            if db == 'blat':
                inputs.append("%s/%s" % (odir, config[db]['x.2bit']))
                inputs.append("%s/%s" % (odir, config[db]['x.ooc']))
            if db == 'gatk':
                inputs.append("%s/%s" % (odir, config[db]['xref']))
                inputs.append("%s/%s" % (odir, config[db]['xref.dict']))
            if db == 'snpeff':
                inputs.append("%s/%s/%s" % (odir, genome, config[db]['xout']))
        if any(x in dbs for x in 'star hisat2 snpeff blastn blastp'.split()):
            inputs.append("%s/%s/15.gff.db" % (genome, config['adir']))
            #inputs.append("%s/%s/25.tandup.tsv" % (genome, config['adir']))
            inputs.append("%s/%s/22.tandup.pro.tsv" % (genome, config['adir']))
            inputs.append("%s/55.rds" % genome)
    return inputs
localrules: all, blast_index, anno1_clean, prepR
rule all:
    input:
        all_inputs

include: "rules/genome.seq.smk"
include: "rules/genome.gff.smk"
rule prepR:
    input:
        chrom_bed = "{genome}/15_intervals/01.chrom.bed",
        chrom_size = "{genome}/15_intervals/01.chrom.sizes",
        gap_bed = "{genome}/15_intervals/11.gap.bed",
        gene_tsv = "{genome}/%s/15.tsv" % config['adir'],
        gene_des = "{genome}/%s/15.desc.tsv" % config['adir'],
    output:
        "{genome}/55.rds"
    params:
        wdir = lambda w: "%s" % w.genome,
    shell:
        """
        genome.prep.R {wildcards.genome}
        """

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
