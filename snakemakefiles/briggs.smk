import os
import os.path as op
from snk.utils import check_config, check_config_ase

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)
check_config_ase(config)

wildcard_constraints:
    sid = "\w+",
    region = ".+(:[0-9]+-[0-9]+)?"
def multiqc_inputs(wildcards):
    inputs = [
        "31_featurecounts/01.txt",
    ]
    t = config['t']
    for i in range(len(t)):
        sid, gt = t['SampleID'][i], t['Genotype'][i]
        inputs.append("21_star/%s/Log.final.out" % sid)
        inputs.append("21_star/%s_unpaired/Log.final.out" % sid)
        inputs.append("logs/trimmomatic/%s.log" % sid)
    return inputs

include: "rules/rename_fq.smk"
if config['paired']:
    include: "rules/trimmomatic_pe.smk"
    include: "rules/star_pe.smk"
else:
    include: "rules/trimmomatic_se.smk"
    include: "rules/star_se.smk"
include: "rules/featurecounts.smk"
include: "rules/bcftools_call_s.smk"
include: "rules/ase.smk"
include: "rules/multiqc.smk"
rule all:
    input:
        "qc/multiqc.html",
        #"25.vcf.gz",
        expand(["22_bam/{sid}.bam"], sid = config['t']['SampleID']),
        #expand(["33_ase/{sid}.tsv"], sid = config['t']['SampleID'])

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))
