import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_reseq2
from jcvi.utils.natsort import natsorted

configfile: 'config.yml'
config = check_config_reseq2(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9]+",
    sid = "[a-zA-Z0-9\-_#]+",
    rid = "[a-zA-Z0-9]+",

localrules: all

def all_outputs(w):
    outputs = []
    for yid, ydic in config['y'].items():
        if not ydic['run']: continue
        odir = "%s/%s" % (config['oid'], yid)
        if config['g'][ydic['ref']]['annotation']:
            outputs.append("%s/%s" % (odir, config['callvnt2']['of38']))
        else:
            outputs.append("%s/%s" % (odir, config['callvnt2']['of37']))
    return outputs

rule all:
    input: all_outputs

include: "rules/callvnt2.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
