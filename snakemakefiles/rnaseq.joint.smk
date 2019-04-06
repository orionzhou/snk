import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_ngs

configfile: 'config.yml'
config = check_config_ngs(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9_]+",
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"

localrules: all, fastq, trimming

def all_outputs(w):
    outputs = []
    for yid in config['y'].keys():
        if not config['y'][yid]['run']: continue
        if config['y'][yid]['meta'] != True:
            outputs.append("%s/data/%s" % (yid, config['merge_trimstats']['out']))
            outputs.append("%s/data/%s" % (yid, config['merge_bamstats']['out']))
        outputs.append("%s/data/%s" % (yid, config['rc2cpm']['out_raw']))
        outputs.append("%s/data/%s" % (yid, config['rc2cpm']['out']))
    return outputs
rule all:
    input: all_outputs

include: "rules/fastq.smk"
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/rnaseq.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
