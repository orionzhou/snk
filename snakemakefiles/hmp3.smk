import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_default

configfile: 'config.yml'
config = check_config_default(config)
workdir: config['dirw']

wildcard_constraints:
    cid = "[0-9]+",
    sid = "[a-zA-Z0-9]+",
    gt = "[a-zA-Z0-9\-_]+",

localrules: all
rule all:
    input:
        config['hmp']['of10'],
        config['hmp']['of21'],
        config['hmp']['of22'],

include: "rules/callvnt.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
