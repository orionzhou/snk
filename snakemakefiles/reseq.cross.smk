import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_ngs_cross

configfile: 'config.yaml'
config = check_config_ngs_cross(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    gt = "[a-zA-Z0-9\-_]+",
    rid = "[a-zA-Z0-9]+",

localrules: all

rule all:
    input:
        "%s/%s" % (config['dird'], config['callvnt_cross']['out']),

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
