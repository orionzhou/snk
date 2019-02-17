import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_rnaseq

configfile: 'config.yaml'
config = check_config_rnaseq(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9]+",
    sid = "[a-zA-Z0-9]+",
    gt = "[a-zA-Z0-9\-_]+",
    rid = "[a-zA-Z0-9]+",

localrules: all, fastq, trimming, merge_trimstats, merge_bamstats

def all_outputs(w):
    outputs = []
    for yid in config['y'].keys():
        pre = "%s/%s" % (yid, config['dird'])
        if config['y'][yid]['meta'] != True:
            outputs.append("%s/%s" % (pre, config['merge_trimstats']['out']))
            outputs.append("%s/%s" % (pre, config['merge_bamstats']['outv']))
            for gt, sids in config['y'][yid]['gt'].items():
                outputs.append("%s/%s/%s.g.vcf.gz" % (pre, config['callvnt']['odir'], gt))
        else:
            outputs.append("%s/%s" % (pre, config['callvnt']['out']))
    return outputs
rule all:
    input: all_outputs

include: "rules/fastq.smk"
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/cleanbam.smk"
include: "rules/callvnt.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
