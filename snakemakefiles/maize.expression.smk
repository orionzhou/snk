import os
import os.path as op
from snk.utils import check_config

configfile: 'config.yaml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"

localrules: all, merge_bamstats

rule all:
    input:
        #expand("%s/{sid}.bam" % config['star']['odir2'], sid = config['SampleID']),
        "%s/%s" % (config['dird'], config['merge_featurecounts']['out']),
        "%s/%s" % (config['dird'], config['merge_bamstats']['out']),
        "%s/%s" % (config['dird'], config['multiqc']['out']),

if config['source'] == 'sra':
    include: "rules/fasterq_dump.smk"
elif config['source'] == 'local_interleaved':
    include: "rules/fq_deinterleave.smk"
elif config['source'] == 'local':
    include: "rules/fq_compress.smk"

if config['readtype'] in ['illumina', 'solid']:
    include: "rules/trimmomatic.smk"
elif config['readtype'] == '3rnaseq':
    include: "rules/bbduk.smk"

if config['mapper'] == 'star':
    include: "rules/star.smk"
elif config['mapper'] == 'hisat2':
    include: "rules/hisat2.smk"

include: "rules/featurecounts.smk"
include: "rules/multiqc.smk"

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
