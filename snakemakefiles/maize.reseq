import os
import os.path as op
from snk.utils import check_config

configfile: 'config.yaml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"
def multiqc_inputs(wildcards):
    inputs = []
    for sid in config['SampleID']:
        inputs.append("%s/%s.txt" % (config['bwa']['odir2'], sid))
        if config['t'][sid]['paired']:
            dirl_trim = 'bbduk_pe' if config['readtype'] == '3rnaseq' else 'trimmomatic_pe'
            inputs.append("%s/%s/%s.log" % (config['dirl'], dirl_trim, sid))
        else:
            dirl_trim = 'bbduk_se' if config['readtype'] == '3rnaseq' else 'trimmomatic_se'
            inputs.append("%s/%s/%s.log" % (config['dirl'], dirl_trim, sid))
    return inputs

localrules: all
if config['source'] == 'sra':
    include: "rules/fasterq_dump.smk"
elif config['source'] == 'local_interleaved':
    include: "rules/fq_deinterleave.smk"
if config['readtype'] in ['illumina', 'solid']:
    include: "rules/trimmomatic.smk"
elif config['readtype'] == '3rnaseq':
    include: "rules/bbduk.smk"
include: "rules/bwa.smk"
include: "rules/multiqc.smk"
rule all:
    input:
        expand("%s/{sid}.bam" % config['bwa']['odir2'], sid = config['SampleID']),
        expand("%s/{sid}.txt" % config['bwa']['odir2'], sid = config['SampleID']),
        "%s/multiqc.html" % config['dird'],

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
