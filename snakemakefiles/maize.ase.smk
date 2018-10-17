import os
import os.path as op
from snk.utils import check_config, check_config_ase

configfile: 'config.yaml'
config = check_config(config)
check_config_ase(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",
    region = ".+(:[0-9]+-[0-9]+)?"

localrules: all, merge_ase

rule all:
    input:
        "%s/%s" % (config['dird'], config['merge_ase']['out']),

include: "rules/ase.smk"

for rule in workflow.rules:
    if rule.name != 'all':
        snakemake.utils.makedirs(op.join(config['dirp'], rule.name))

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
