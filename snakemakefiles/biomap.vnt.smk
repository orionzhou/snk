import os
import os.path as op
from snk.utils import get_resource
from snk.utils import check_config_default

configfile: 'config.yml'
config = check_config_default(config)
workdir: config['dirw']

wildcard_constraints:
    sid = "[a-zA-Z0-9]+",

localrules: all
rule all:
    input:
        config['biomap']['of04'],
        config['biomap']['of05'],

rule liftover:
    input:
        config['biomap']['if01'],
    output:
        vcf = temp(config['biomap']['of03']),
        tbi = temp(config['biomap']['of03'] + ".tbi"),
    params:
        refn = config['biomap']['refn'],
        chain4 = config['biomap']['chain4'],
        rej = config['biomap']['of03'].replace('.vcf.gz','.rej.vcf'),
        N = "%s" % config['biomap']['liftover']['id'],
        e = "%s/%s.e" % (config['dirp'], config['biomap']['liftover']['id']),
        o = "%s/%s.o" % (config['dirp'], config['biomap']['liftover']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'biomap','liftover')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'biomap','liftover')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'biomap','liftover')['mem']
    threads: config['biomap']['liftover']['ppn']
    conda: "envs/gatk.yml"
    shell:
        """
        gatk --java-options "-Xmx{params.mem}G" LiftoverVcf \
        --MAX_RECORDS_IN_RAM 5000 \
        --TMP_DIR /scratch.global/zhoux379/temp/ \
        --USE_JDK_DEFLATER --USE_JDK_INFLATER \
        -I {input} -O {output.vcf} --REJECT {params.rej} \
        -C {params.chain4} -R {params.refn}
        """

rule filter:
    input: config['biomap']['of03']
    output:
        snp = config['biomap']['of04'],
        idl = config['biomap']['of05'],
    params:
        tmp = config['tmpdir'],
        snp_stat = lambda w, output: output.snp.replace(".vcf.gz", ".stats.txt"),
        idl_stat = lambda w, output: output.idl.replace(".vcf.gz", ".stats.txt"),
        snp_sites = lambda w, output: output.snp.replace(".vcf.gz", ".sites.vcf.gz"),
        idl_sites = lambda w, output: output.idl.replace(".vcf.gz", ".sites.vcf.gz"),
        snp_bed = lambda w, output: output.snp.replace(".vcf.gz", ".bed"),
        idl_bed = lambda w, output: output.idl.replace(".vcf.gz", ".bed"),
        N = config['biomap']['filter']['id'],
        e = "%s/%s.e" % (config['dirp'], config['biomap']['filter']['id']),
        o = "%s/%s.o" % (config['dirp'], config['biomap']['filter']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'biomap','filter')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'biomap','filter')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'biomap','filter')['mem']
    threads: config['biomap']['filter']['ppn']
    conda: "envs/samtools.yml"
    shell:
#        bcftools view -m2 -M2 -v indels -i 'MAF>0.05' {input} -Oz -o {output.idl}
        """
        bcftools view -m2 -M2 -v snps -i 'MAF>0.05' {input} -Oz -o {output.snp}
        bcftools view -m2 -c 4 -v indels -i 'MAF>0.05' {input} -Oz -o {output.idl}
        bcftools index -t {output.snp}
        bcftools index -t {output.idl}
        bcftools stats -s - {output.snp} > {params.snp_stat}
        bcftools stats -s - {output.idl} > {params.idl_stat}
        bcftools view -G {output.snp} -Ou | bcftools annotate -x INFO -Oz -o {params.snp_sites}
        bcftools view -G {output.idl} -Ou | bcftools annotate -x INFO -Oz -o {params.idl_sites}
        bcftools index -t {params.snp_sites}
        bcftools index -t {params.idl_sites}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.snp} -o {params.snp_bed}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.idl} -o {params.idl_bed}
        """


onsuccess:
    shell("mail -s 'Success: %s' %s < {log}" % (config['dirw'], config['email']))
onerror:
    shell("mail -s 'Error: %s' %s < {log}" % (config['dirw'], config['email']))
