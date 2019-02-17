import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_default

configfile: 'config.yaml'
config = check_config_default(config)
workdir: config['dirw']

wildcard_constraints:
    cid = "[0-9]+",
    sid = "[a-zA-Z0-9]+",
    gt = "[a-zA-Z0-9\-_]+",

localrules: all
rule all:
    input:
        config['hmp']['merge']['out']

rule liftover:
    input: "%s/merged_flt_c{cid}.imputed.vcf.gz" % config['hmp']['liftover']['idir']
    output:
        vcf = "%s/{cid}.vcf.gz" % config['hmp']['liftover']['odir'],
        tbi = "%s/{cid}.vcf.gz.tbi" % config['hmp']['liftover']['odir'],
    params:
        refo = config['hmp']['refo'],
        refn = config['hmp']['refn'],
        chain3 = config['hmp']['chain3'],
        chain4 = config['hmp']['chain4'],
        o1 = "%s/{cid}.1.vcf" % config['hmp']['liftover']['odir'],
        o2 = "%s/{cid}.2.vcf.gz" % config['hmp']['liftover']['odir'],
        o3 = "%s{cid}.3.vcf.gz" % config['hmp']['liftover']['odir'],
        o4 = "%s/{cid}.4.vcf.gz" % config['hmp']['liftover']['odir'],
        o5 = "%s/{cid}.5.vcf.gz" % config['hmp']['liftover']['odir'],
        rej = "%s/{cid}.rej.vcf" % config['hmp']['liftover']['odir'],
        N = "%s.{cid}" % config['hmp']['liftover']['id'],
        e = "%s/%s/{cid}.e" % (config['dirp'], config['hmp']['liftover']['id']),
        o = "%s/%s/{cid}.o" % (config['dirp'], config['hmp']['liftover']['id']),
        mem = lambda w, resources: resources.mem - 2
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','liftover')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','liftover')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','liftover')['mem']
    threads: config['hmp']['liftover']['ppn']
    shell:
        """
        set +u
        source $soft/miniconda3/etc/profile.d/conda.sh

        conda activate crossmap
        CrossMap.py vcf {params.chain3} {input} {params.refo} {params.o1}

        conda activate samtools
        bcftools sort -m {params.mem}G -Oz -o {params.o2} {params.o1}
        bcftools norm -c x -m -any -f {params.refo} {params.o2} -o {params.o3}
        bcftools filter -e "REF==ALT[0]" {params.o3} -o {params.o4}
        bcftools norm -m +any -f {params.refo} {params.o4} {params.o5}

        conda activate gatk
        picard -Xmx{params.mem}G LiftoverVcf \
        I={params.o5} O={output.vcf} REJECT={params.rej} \
        C={params.chain4} R={params.refn}
        """

rule merge:
    input: expand("%s/{cid}.vcf.gz" % config['hmp']['liftover']['odir'], cid=range(1,11))
    output: config['hmp']['merge']['out']
    params:
        o1 = "10.unsorted.vcf.gz",
        N = config['hmp']['merge']['id'],
        e = "%s/%s.e" % (config['dirp'], config['hmp']['merge']['id']),
        o = "%s/%s.o" % (config['dirp'], config['hmp']['merge']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','merge')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','merge')['mem']
    threads: config['hmp']['merge']['ppn']
    shell:
        """
        set +u
        source $soft/miniconda3/etc/profile.d/conda.sh

        conda activate samtools
        bcftools concat -a -n -Oz -o {params.o1} {input}
        bcftools sort -m {params.mem}G -Oz -o {output} {params.o1}
        """

onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
