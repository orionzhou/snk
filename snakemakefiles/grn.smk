import os
import os.path as op
import pandas as pd
import yaml
from snakemake.utils import update_config, makedirs
from snk.utils import get_resource, check_config_default

def check_config(c):
    c = check_config_default(c)
    c['dirw'] = c['dirc']

    f_cfg = op.join(c['dird'], c['grn']['cfg'])
    df = pd.read_excel(f_cfg, sheet_name=0, header=0)
    c['grn']['evtype'] = c['grn']['evtype'].split()
    c['nid'] = []
    c['t'] = dict()
    for i in range(len(df)):
        nid = df['nid'][i]
        if not nid.startswith('np_gibeerish'):
            c['nid'].append(nid)
            c['t'][nid] = {x: df[x][i] for x in list(df)}
    return c

configfile: 'config.yml'
config = check_config(config)
workdir: config['dirw']

wildcard_constraints:
    nid = "[a-zA-Z0-9_]+",
    a_nid = "[a-zA-Z0-9_]+",
    b_nid = "[a-zA-Z0-9_]+",
    evtype = "[a-zA-Z0-9]+",

rule all:
    input:
        expand("%s/{nid}.tsv" % config['grn']['od11'], nid=config['nid']),
        expand("%s/{nid}.pkl" % config['grn']['od14'], nid=config['nid']),
        expand("%s/{nid}.rds" % config['grn']['od14'], nid=config['nid']),
        expand("%s/01.meval.rds" % config['dirr'], nid=config['nid']),
#        expand("%s/01.{evtype}.rds" % config['dirr'], evtype=config['grn']['evtype'])

include: "rules/grn.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
