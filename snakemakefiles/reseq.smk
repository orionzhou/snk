import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_ngs

configfile: 'config.yml'
config = check_config_ngs(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9]+",
    sid = "[a-zA-Z0-9\-_]+",
    gt = "[a-zA-Z0-9\-_]+",
    rid = "[a-zA-Z0-9]+",

localrules: all, fastq, trimming, merge_trimstats, merge_bamstats

def all_outputs(w):
    outputs = []
    for yid in config['y'].keys():
        if not config['y'][yid]['run']: continue
        if config['y'][yid]['meta'] != True:
            outputs.append("%s/%s" % (yid, config['cleanbam']['of24c']))
            outputs.append("%s/data/%s" % (yid, config['merge_trimstats']['out']))
            outputs.append("%s/data/%s" % (yid, config['merge_bamstats']['outv']))
#            for gt, sids in config['y'][yid]['gt'].items():
#                outputs.append("%s/%s/%s.g.vcf.gz" % (config['callvnt']['odir'], yid, gt))
        else:
            if config['g'][config['y'][yid]['reference']]['annotation']:
                outputs.append("%s/%s" % (yid, config['callvnt']['of38']))
            else:
                outputs.append("%s/%s" % (yid, config['callvnt']['of37']))
    return outputs
rule all:
    input: all_outputs

include: "rules/fastq.smk"
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/bcftools.smk"
include: "rules/cleanbam.smk"
include: "rules/callvnt.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
