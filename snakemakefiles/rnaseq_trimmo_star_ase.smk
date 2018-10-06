import os
import os.path as op
from snk.utils import check_config, check_config_ase

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)
#check_config_ase(config)

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"
def multiqc_inputs(wildcards):
    inputs = []
    for sid in config['SampleID']:
        inputs.append("%s/%s.txt.summary" % (config['featurecounts']['odir'], sid))
        if config['t'][sid]['paired']:
            inputs.append("%s/%s_p/Log.final.out" % (config['star']['odir1'], sid))
            inputs.append("%s/%s_u/Log.final.out" % (config['star']['odir1'], sid))
            inputs.append("logs/trimmomatic_pe/%s.log" % sid)
        else:
            inputs.append("%s/%s/Log.final.out" % (config['star']['odir1'], sid))
            inputs.append("logs/trimmomatic_se/%s.log" % sid)
    return inputs

#include: "rules/fq_rename.smk"
include: "rules/fq_deinterleave.smk"
include: "rules/trimmomatic.smk"
include: "rules/star.smk"
include: "rules/featurecounts.smk"
include: "rules/multiqc.smk"
rule all:
    input:
        "%s/multiqc.html" % config['dird'],
        "%s/%s" % (config['dird'], config['featurecounts']['out']),
        expand("%s/{sid}.bam" % config['star']['odir2'], sid = config['SampleID']),
        #expand(["33_ase/{sid}.tsv"], sid = config['SampleID'])

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
