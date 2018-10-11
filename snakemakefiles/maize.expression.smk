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
        inputs.append("%s/%s.txt.summary" % (config['featurecounts']['odir'], sid))
        if config['mapper'] == 'hisat2':
            inputs.append("%s/%s.txt" % (config['hisat2']['odir1'], sid))
        if config['t'][sid]['paired']:
            if config['mapper'] == 'star':
                inputs.append("%s/%s_p/Log.final.out" % (config['star']['odir1'], sid))
                inputs.append("%s/%s_u/Log.final.out" % (config['star']['odir1'], sid))
            dirl_trim = 'bbduk_pe' if config['readtype'] == '3rnaseq' else 'trimmomatic_pe'
            inputs.append("%s/%s/%s.log" % (config['dirl'], dirl_trim, sid))
        else:
            if config['mapper'] == 'star':
                inputs.append("%s/%s/Log.final.out" % (config['star']['odir1'], sid))
            dirl_trim = 'bbduk_se' if config['readtype'] == '3rnaseq' else 'trimmomatic_se'
            inputs.append("%s/%s/%s.log" % (config['dirl'], dirl_trim, sid))
    return inputs

localrules: all
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
rule all:
    input:
        #expand("%s/{sid}.bam" % config['star']['odir2'], sid = config['SampleID']),
        "%s/%s" % (config['dird'], config['featurecounts']['out']),
        "%s/multiqc.html" % config['dird'],

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
