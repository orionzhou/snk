import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_reseq
from jcvi.utils.natsort import natsorted
import pandas as pd

configfile: 'config.yml'
config = check_config_reseq(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9]+",
    sid = "[a-zA-Z0-9\-]+",
    gt = "[a-zA-Z0-9\-_]+",
    rid = "[a-zA-Z0-9]+",

localrules: all, trimming, merge_trimstats, bam9_merge_stats, cv14_vcfstats

def all_outputs(w):
    outputs = []
    for yid, ydic in config['y'].items():
        if not ydic['runD']: continue
        odir = "%s/%s" % (config['oid'], yid)
        gts = ydic['Genotypes']
        sufs = '.g.vcf.gz .g.vcf.gz.tbi .txt'.split()
#        outputs.append("%s/%s" % (odir, config['trimming']['out']))
#        outputs.append("%s/%s" % (odir, config['merge_bamstats']['outv']))
        outputs += expand("%s/%s/{gt}{suf}" % (config['callvnt']['odir'], yid), gt = gts, suf = sufs)
        outputs.append("%s/%s" % (odir, config['callvnt']['out']))
    return outputs

rule all:
    input: all_outputs

include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/bcftools.smk"
include: "rules/cleanbam.smk"
include: "rules/callvnt1.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
