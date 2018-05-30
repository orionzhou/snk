def gatk_hc_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['haplotype_caller']:
        java_mem = config['gatk']['haplotype_caller']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
   
    java_option = "-Xmx%s -Djava.io.tmpdir=%s" % (java_mem, java_tmpdir)
    cmd = "%s --java-options \"%s\" %s" % (config['gatk']['cmd'], 
            java_option, config['gatk']['haplotype_caller']['mode'])
    return cmd

rule gatk_haplotype_caller:
    input:
        expand(["%s/{sid}_{gt}.bam" % config['gatk']['idir']],
            zip, sid = config['t']['sid'], gt = config['t']['genotype'])
    output:
        "%s/01.vcf" % config['gatk']['odir']
    params:
        cmd = gatk_hc_cmd,
        ref = config['gatk']['ref'],
        extra = config['gatk']['haplotype_caller']['extra'],
    threads:
        config["gatk"]['haplotype_caller']["threads"]
    run:
        shell("{params.cmd} "
            "-R {params.ref} "
            "{params.extra} "
            "-I {input} "
            "-O {output}")
