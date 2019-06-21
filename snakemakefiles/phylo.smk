import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_phylo

configfile: 'config.yml'
config = check_config_phylo(config)
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
        outputs.append("%s/%s" % (odir, config['phylo']['of35']))
    return outputs
rule all:
    input: all_outputs

include: "rules/phylo.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
