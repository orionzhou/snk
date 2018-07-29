def gatk_hc_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['haplotype_caller']:
        java_mem = config['gatk']['haplotype_caller']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
    cmd = "%s -Xmx%s -Djava.io.tmpdir=%s" % (config['gatk']['cmd'], java_mem, java_tmpdir)
    return cmd

rule gatk_haplotype_caller:
    input:
        "%s/{sid}.bam" % config['gatk']['idir']
    output:
        temp("%s/{sid}/{region}.g.vcf.gz" % config['gatk']['odir'])
    log:
        "%s/gatk/{sid}/{region}.log" % config['dirl']
    params:
        cmd = gatk_hc_cmd,
        ref = config['gatk']['ref'],
        extra = config['gatk']['haplotype_caller']['extra'],
    threads:
        config["gatk"]['haplotype_caller']["threads"]
    shell:
        #-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation 
        """
        source activate gatk

        {params.cmd} -T HaplotypeCaller \
        -nct {threads} \
        -R {params.ref} \
        -ERC GVCF \
        -G Standard -G AS_Standard \
        -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
        -U ALLOW_N_CIGAR_READS \
        -L {wildcards.region} \
        {params.extra} \
        -I {input} \
        -o {output} 2>{log}
        """

def bcftools_inputs(wildcards):
    sid = wildcards.sid
    inputs = ["%s/%s/%s.g.vcf.gz" % 
        (config['gatk']['odir'], sid, x) 
        for x in config['gatk']['bcftools_concat']['regions'].split(" ")]
    return inputs

rule bcftools_concat:
    input:
        bcftools_inputs
    output:
        temp("%s/{sid}.g.vcf.gz" % config['gatk']['odir'])
    log:
        "%s/gatk/{sid}.log" % config['dirl']
    params:
        extra = config['gatk']['bcftools_concat']['extra']
    threads:
        config["gatk"]['bcftools_concat']["threads"]
    shell:
        "bcftools concat --threads {threads} {input} -O z -o {output}"
