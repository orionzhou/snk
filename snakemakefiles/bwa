import os
import os.path as op
from snk.utils import check_config

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"
def multiqc_inputs(wildcards):
    inputs = []
    for sid in config['SampleID']:
        if config['t'][sid]['paired']:
            inputs.append("%s/trimmomatic_pe/%s.log" % (config['dirl'], sid))
        else:
            inputs.append("%s/trimmomatic_se/%s.log" % (config['dirl'], sid))
        inputs.append("%s/%s.txt" % (config['bwa']['odir2'], sid))
    return inputs

#include: "rules/fasterq_dump.smk"
include: "rules/fq_deinterleave.smk"
include: "rules/trimmomatic.smk" 
include: "rules/bwa.smk" 
include: "rules/callvnt_bcftools.smk" 
include: "rules/multiqc.smk"
rule all:
    input:
        "%s/multiqc.html" % config['dird'],
        config['callvnt']['outfile'],
        expand("%s/{sid}.g.vcf.gz" % config['callvnt']['odir1'], sid = config['SampleID'])

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
