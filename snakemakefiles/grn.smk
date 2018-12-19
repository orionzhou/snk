import os
import os.path as op
import pandas as pd
import yaml
from snakemake.utils import update_config, makedirs
from snk.utils import get_resource

def check_config(c):
    fy = open(c['config_default'], 'r')
    config_default = yaml.load(fy)
    update_config(config_default, c)
    c = config_default
    
    for subdir in [c['dirw'], c['tmpdir']]:
        if not op.isdir(subdir):
            makedirs(subdir)
    
    for rsubdir in [c['dirp']]:
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)
    
    df = pd.read_excel(c['grn']['cfg'], sheet_name=0, header=0)
    c['grn']['evtype'] = c['grn']['evtype'].split()
    c['nid'] = []
    c['t'] = dict()
    for i in range(len(df)):
        nid = df['nid'][i]
        if not nid.startswith('np_gibeerish'):
            c['nid'].append(nid)
            c['t'][nid] = {x: df[x][i] for x in list(df)}
    return c

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)

wildcard_constraints:
    study = "[a-zA-Z0-9]+",
    evtype = "[a-zA-Z0-9]+",

rule all:
    input:
        expand("%s/{nid}.rda" % config['grn']['pkl2rda']['odir'], nid=config['nid']),
        expand("%s/01.{evtype}.rds" % config['grn']['eval_merge']['odir'], evtype=config['grn']['evtype'])

include: "rules/grn.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
