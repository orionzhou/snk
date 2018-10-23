import os
import os.path as op
from snakemake.utils import update_config, makedirs
from astropy.table import Table, Column
import yaml
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
    
    t = Table.read(c['genie3']['cfg'], format = 'ascii.tab')
    c['nid'] = t['nid']
    c['t'] = dict()
    cols = t.colnames
    for i in range(len(t)):
        nid = t['nid'][i]
        if nid not in ['n16b','n16c','n99b_1','nc01','nc02','nc03']:
            continue
        c['t'][nid] = {x: t[x][i] for x in cols}
    
    return c

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)

wildcard_constraints:
    study = "[a-zA-Z0-9]+",

rule all:
    input:
        expand(["%s/{nid}.pkl" % config['genie3']['odir']], nid = config['nid']),
        #expand(["%s/{nid}.pkl" % config['genie3']['odir']], nid = ['n13a','n13b']),

include: "rules/genie3.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
