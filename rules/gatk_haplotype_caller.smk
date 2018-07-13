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
        "%s/{sample}.bam" % config['gatk']['idir']
    output:
        "%s/{sample}.g.vcf.gz" % config['gatk']['odir']
    log:
        "%s/gatk/{sample}.log" % config['dirl']
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
        {params.extra} \
        -I {input} \
        -o {output} 2>{log}
        """
