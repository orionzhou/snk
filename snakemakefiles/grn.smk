import os
import os.path as op
import pandas as pd
from snk.utils import get_resource
from snk.utils import check_config_grn

configfile: 'config.yml'
config = check_config_grn(config)
workdir: config['dirw']

wildcard_constraints:
    nid = "[a-zA-Z0-9_]+",
    a_nid = "[a-zA-Z0-9_]+",
    b_nid = "[a-zA-Z0-9_]+",
    gopt = "[a-zA-Z0-9]+",
    eopt = "[a-zA-Z0-9]+",

rule all:
    input:
        expand("%s/{nid}.tsv" % config['grn']['od12'], nid=config['nid']),
        expand("%s/{gopt}.{nid}.rds" % config['grn']['od15'], gopt=config['grn']['gopts'], nid=config['nid']),
        expand("%s/{gopt}.rds" % config['oid'], gopt = config['grn']['gopts']),
        expand("%s/{gopt}.{eopt}.rds" % config['oid'], gopt=config['grn']['gopts'], eopt=config['grn']['eopts'])

include: "rules/grn.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
