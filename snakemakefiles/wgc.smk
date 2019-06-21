import os
import os.path as op
import pandas as pd
import yaml
from snk.utils import get_resource
from snk.utils import check_config_wgc

configfile: 'config.yml'
config = check_config_wgc(config)
workdir: config['dirw']

wildcard_constraints:
    qry = "[a-zA-Z0-9]+",
    tgt = "[a-zA-Z0-9]+",
    idx = "[0-9]+",
    tchrom = "[a-zA-Z0-9]+"

localrules: all

def all_inputs(wildcards):
    inputs = []
    inputs.append("%s/%s" % (config['wgc']['od50'], config['wgc']['of50']))
    return inputs
    for qry, tgt in config['comps']:
        odir = "%s/%s_%s" % (config['dirr'], qry, tgt)
        inputs.append("%s/%s" % (odir, config['wgc']['rf10']))
        inputs.append("%s/15.%s.tsv" % (odir, tgt))
        if config['x'][qry]['annotation']:
            inputs.append("%s/15.%s.tsv" % (odir, qry))
            inputs.append("%s/%s/%s" % (odir, config['wgc']['rd20'], config['wgc']['out']))
    return inputs
rule all:
    input:
        all_inputs

include: "rules/wgc.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
