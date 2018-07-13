def gatk_gg_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['genotype_gvcfs']:
        java_mem = config['gatk']['genotype_gvcfs']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
    cmd = "%s -Xmx%s -Djava.io.tmpdir=%s" % (config['gatk']['cmd'], java_mem, java_tmpdir)
    return cmd

rule gatk_genotype_gvcfs:
    input:
        expand(["%s/{sid}_{gt}.g.vcf.gz" % config['gatk']['odir']],
            zip, sid = config['t']['SampleID'], gt = config['t']['Genotype'])
    output:
        vcf = config['gatk']['outfile']
    log:
        "%s/gatk.log" % config['dirl']
    params:
        cmd = gatk_gg_cmd,
        ref = config['gatk']['ref'],
        gvcfs = expand(["-V %s/{sid}_{gt}.g.vcf.gz" % config['gatk']['odir']],
            zip, sid = config['t']['SampleID'], gt = config['t']['Genotype']),
        gvcf = "%s/all.g.vcf.gz" % config['gatk']['odir'],
        extra = config['gatk']['genotype_gvcfs']['extra'],
    threads:
        config["gatk"]['genotype_gvcfs']["threads"]
    shell:
        """
        source activate gatk

        {params.cmd} -T CombineGVCFs \
        -R {params.ref} \
        {params.extra} \
        {params.gvcfs} \
        -o {params.gvcf} 2>{log}

        {params.cmd} -T GenotypeGVCFs \
        -nt {threads} \
        -R {params.ref} \
        {params.extra} \
        -V {params.gvcf} \
        -o {output.vcf} 2>{log}
        """
