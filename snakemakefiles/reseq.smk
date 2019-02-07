import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_ngs

configfile: 'config.yaml'
config = check_config_ngs(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    gt = "[a-zA-Z0-9\-_]+",
    rid = "[a-zA-Z0-9]+",

localrules: all, merge_bamstats, merge_trimstats

rule all:
    input:
        expand("%s/%s/{gt}.g.vcf.gz" % (config['dird'], config['callvnt']['odir1']), gt = config['Genotypes']),
#        "%s/%s" % (config['dird'], config['merge_trimstats']['out']),
#        "%s/%s" % (config['dird'], config['merge_bamstats']['out']),
#        "%s/%s" % (config['dird'], config['callvnt']['out']),

if config['source'] == 'sra':
    include: "rules/fasterq_dump.smk"
elif config['source'] == 'local_interleaved':
    include: "rules/fq_deinterleave.smk"
elif config['source'] == 'local':
    include: "rules/fq_compress.smk"

include: "rules/fastp.smk"
include: "rules/bwa.smk"
include: "rules/cleanbam.smk"
include: "rules/callvnt_gatk.smk"
include: "rules/report.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
