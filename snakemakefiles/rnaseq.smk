import os
import os.path as op
from snk.utils import check_config_ngs
from snk.utils import get_resource

configfile: 'config.yaml'
config = check_config_ngs(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"

localrules: all, merge_bamstats, merge_trimstats

rule all:
    input:
        "%s/%s" % (config['dird'], config['merge_featurecounts']['out']),
        "%s/%s" % (config['dird'], config['merge_trimstats']['out']),
        "%s/%s" % (config['dird'], config['merge_bamstats']['out']),
        "%s/%s" % (config['dird'], config['multiqc']['out']),

if config['source'] == 'sra':
    include: "rules/fasterq_dump.smk"
elif config['source'] == 'local_interleaved':
    include: "rules/fq_deinterleave.smk"
elif config['source'] == 'local':
    include: "rules/fq_compress.smk"

include: "rules/fastp.smk"
if config['mapper'] == 'star':
    include: "rules/star.smk"
elif config['mapper'] == 'hisat2':
    include: "rules/hisat2.smk"

include: "rules/featurecounts.smk"
include: "rules/report.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
