import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_brbseq

configfile: 'config.yml'
config = check_config_brbseq(config)
workdir: config['dirw']

wildcard_constraints:
    yid = "[a-zA-Z0-9_]+",
    sid = "[a-zA-Z0-9]+",

localrules: all

sids = ['batch1','batch2a','batch2b','batch3']
def all_outputs(w):
    outputs = []
    for sid in sids:
        odir = op.join(config['oid'], sid)
        outputs.append("%s/%s/output.dge.umis.txt" % (config['brbseq']['od21'], sid))
#        outputs.append("%s/%s/A1.tsv" % (config['brbseq']['od04'], sid))
    return outputs
rule all:
    input: all_outputs

include: "rules/brbseq.smk"

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
