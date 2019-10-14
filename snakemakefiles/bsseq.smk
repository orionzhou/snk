import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_bsseq

configfile: 'config.yml'
config = check_config_bsseq(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    mid = "[a-zA-Z0-9]+",
    part = "[a-z]+",
    gt = "[a-zA-Z0-9\-_]+",
    rid = "[a-zA-Z0-9]+",

localrules: all, trimming, merge_trimstats, merge_bamstats

def all_outputs(w):
    outputs = []
    for yid, ydic in config['y'].items():
        if not ydic['runB']: continue
        sids = config['y'][yid]['SampleID']
#        outputs += expand("%s/%s/{sid}.rds" % (config['oid'], yid), sid=sids)
        mids = config['y'][yid]['MergeID']
        outputs += expand("%s/%s/{mid}.cx.gz" % (config['oid'], yid), mid=mids)
        outputs += expand("%s/%s/{mid}.bed.gz" % (config['oid'], yid), mid=mids)
    return outputs

rule all:
    input: all_outputs

include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/cleanbam.smk"
include: "rules/bsseq.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
