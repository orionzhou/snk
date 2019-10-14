import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_reseq3
from jcvi.utils.natsort import natsorted

configfile: 'config.yml'
config = check_config_reseq3(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    yid = "[a-zA-Z0-9]+",

localrules: all
def all_outputs(w):
    outputs = []
    for yid, ydic in config['y'].items():
        if not ydic['run']: continue
        odir = "%s/%s" % (config['oid'], yid)
        outputs.append("%s/%s" % (yid, config['callvnt3']['of11']))
        outputs += expand("%s/%s/{gt}.bcf" % (yid, config['callvnt3']['od22']), gt = ydic['ase_genotypes'])
    return outputs
rule all:
    input: all_outputs

include: "rules/callvnt3.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
