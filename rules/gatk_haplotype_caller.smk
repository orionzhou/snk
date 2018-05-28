def gatk_hc_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['haplotype_caller']:
        java_mem = config['gatk']['haplotype_caller']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
    
    cmd = "%s -T %s -Xmx%s -Djava.io.tmpdir=%s" % (config['gatk']['cmd'], 
        config['gatk']['haplotype_caller']['mode'],
        java_mem, java_tmpdir)
    return cmd

rule gatk_haplotype_caller:
    input:
        expand(["22.bam/{sid}_{gt}.bam"],
            zip, sid = config['t']['sid'], gt = config['t']['genotype'])
    output:
        "25.gatk/01.vcf"
    params:
        cmd = gatk_hc_cmd,
        ref = config['gatk']['ref'],
        extra = config['gatk']['haplotype_caller']['extra'],
    threads:
        config["gatk"]['haplotype_caller']["threads"]
    run:
        shell("{params.cmd} "
        "-nct {threads} "
	"-R {params.ref} "
	"-I {input} "
        "{params.extra} "
        "-o {output} "
        "{input}")

