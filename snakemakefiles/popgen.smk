import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_popgen

configfile: 'config.yml'
config = check_config_popgen(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9_-]+",
    gid = "[a-zA-Z0-9_]+",

localrules: all
def all_outputs(w):
    outputs = []
    outputs.append(config['popgen']['of13'])
    #outputs += expand("%s/{sid}.fas" % config['popgen']['od06'], sid = config['sids'])
    return outputs
rule all:
    input: all_outputs

include: "rules/popgen.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
