import os
import os.path as op
from snk.utils import get_resource, check_config_barn

configfile: 'config.yml'
config = check_config_barn(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9]+",
    sid = "[a-zA-Z0-9\-]+",

localrules: all
def all_outputs(w):
    outputs = []
    for yid, ydic in config['y'].items():
        if not ydic['run']: continue
        odir = "%s/%s" % (config['barn']['fqdir'], yid)
        fqs = expand("%s/{sid}{suf}" % odir, sid=ydic['SampleID'], suf=['.fq.gz', '_1.fq.gz', '_2.fq.gz'])
        outputs.append("%s/%s/%s.tsv" % (config['dirh'], config['barn']['odir'], yid))
    return outputs
rule all:
    input: all_outputs

include: "rules/barn.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
