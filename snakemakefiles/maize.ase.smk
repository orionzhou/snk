import os
import os.path as op
from snk.utils import check_config_ngs
from snk.utils import get_resource

def check_config_ase(c):
    for subdir in [c['ase']['vdir']]: 
        if not op.isdir(subdir):
            mkdir(subdir)
    t = c['t']
    c['vcf'] = dict()
    c['vbed'] = dict()
    for sid in c['SampleID']:
        gt = t[sid]['Genotype']
        fv = op.join(c['ase']['vdir'], "%s.vcf" % gt)
        fb = op.join(c['ase']['vdir'], "%s.bed" % gt)
        if not op.isfile(fb):
            fv = op.join(c['ase']['vdir2'], "%s.vcf" % gt)
            fb = op.join(c['ase']['vdir2'], "%s.bed" % gt)
        #assert op.isfile(fv), "no vcf found: %s" % fv
        assert op.isfile(fb), "no variant-bed found: %s" % fb
        c['vcf'][sid] = fv
        c['vbed'][sid] = fb

configfile: 'config.yaml'
config = check_config_ngs(config)
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

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
