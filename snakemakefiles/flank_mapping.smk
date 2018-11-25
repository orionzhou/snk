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
    
    for rsubdir in [c['dirl'], c['dirp']]:
        subdir = op.join(c['dirw'], rsubdir)
        if not op.isdir(subdir):
            makedirs(subdir)

    t = Table.read(c['cfg'], format = 'ascii.tab')
    c['t'] = dict()
    c['pairs'] = set()
    for i in range(len(t)):
        sid = t['sid'][i]
        type = t['type'][i]
        genotype = t['genotype'][i]
        tgt = t['tgt'][i]
        opt	= t['opt'][i]
        mode = t['mode'][i]
        c['t'][sid] = {
                'type': type,
                'genotype': genotype,
                'tgt': tgt,
                'opt': opt,
                'mode': mode
        }
        c['pairs'].add((genotype, tgt))
    return c

configfile: 'config.yaml'
workdir: config['dirw']
config = check_config(config)

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    type = "[a-zA-Z0-9]+",
    genotype = "[a-zA-Z0-9]+",
    tgt = "[a-zA-Z0-9]+",
    mode = "[a-zA-Z0-9]+",
    opt = "[a-zA-Z0-9\+]+",

localrules: all, fm1_get_bed, fm2_get_seq, fm2_get_seq_se, fm2_get_seq_merged
rule all:
    input:
        expand("%s/{sid}.filtered.tsv" % config['fm']['odir4'], sid = config['t'].keys()),
        #expand("%s/{p[0]}_{p[1]}.tsv" % config['lastz']['odir'], p = config['pairs'])
        "%s/PH207_B73.tsv" % config['lastz']['odir'],
        "%s/PH207_W22.tsv" % config['lastz']['odir'],
        "%s/PH207_PHB47.tsv" % config['lastz']['odir'],
        "%s/PH207_Mo17.tsv" % config['lastz']['odir'],
        "%s/W22_B73.tsv" % config['lastz']['odir'],
        "%s/W22_PH207.tsv" % config['lastz']['odir'],
        "%s/W22_Mo17.tsv" % config['lastz']['odir'],
        "%s/W22_PHB47.tsv" % config['lastz']['odir'],

include: "rules/polyte_fm.smk"
include: "rules/polyte_lastz.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
