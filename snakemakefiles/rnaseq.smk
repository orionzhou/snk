import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_rnaseq

def normalize_genotype(gt):
    gts = gt.split('x')
    if len(gts) == 2 and gts[0] > gts[1]:
        return 'x'.join((gts[1], gts[0]))
    else:
        return gt

configfile: 'config.yml'
config = check_config_rnaseq(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9_]+",
    sid = "[a-zA-Z0-9]+",
    gt = "[a-zA-Z0-9\-_]+",
    region = ".+(:[0-9]+-[0-9]+)?"

localrules: all, trimming, merge_trimstats, merge_bamstats

def all_outputs(w):
    outputs = []
    for yid in config['y'].keys():
        if not config['y'][yid]['runR']: continue
        odir = op.join(config['oid'], yid)
        outputs.append("%s/%s" % (odir, config['trimming']['out']))
        outputs.append("%s/%s" % (odir, config['merge_bamstats']['out']))
        outputs.append("%s/%s" % (odir, config['rnaseq']['out_fcnt']))
#        outputs.append("%s/%s" % (odir, config['rnaseq']['out_mmq']))
#        if config['y'][yid]['ase']:
#            outputs.append("%s/%s" % (odir, config['rnaseq']['out_ase']))
        outputs.append("%s/%s" % (odir, config['rnaseq']['out_rcpm']))
        outputs.append("%s/%s" % (odir, config['rnaseq']['out_cpm']))
    return outputs
rule all:
    input: all_outputs

include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/rnaseq.smk"
include: "rules/ase.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
