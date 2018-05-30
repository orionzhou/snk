def gatk_vf_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['variant_filtration']:
        java_mem = config['gatk']['variant_filtration']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
   
    java_option = "-Xmx%s -Djava.io.tmpdir=%s" % (java_mem, java_tmpdir)
    cmd = "%s --java-options \"%s\"" % (config['gatk']['cmd'], java_option)
    return cmd

rule gatk_variant_filtration:
    input:
        "%s/01.vcf" % config['gatk']['odir']
    output:
        "%s/10.vcf" % config['gatk']['odir']
    params:
        pre = config['gatk']['odir'],
        cmd = gatk_vf_cmd,
        ref = config['gatk']['ref'],
        extra = config['gatk']['variant_filtration']['extra'],
    run:
        shell("{params.cmd} SelectVariants -R {params.ref} "
            "-V {params.pre}/01.vcf "
            "--select-type-to-include SNP "
            "-O {params.pre}/02.snp.vcf"),
        shell("{params.cmd} SelectVariants -R {params.ref} "
            "-V {params.pre}/01.vcf "
            "--select-type-to-include INDEL "
            "-O {params.pre}/02.indel.vcf"),
        shell("{params.cmd} VariantFiltration -R {params.ref} "
            "-V {params.pre}/02.snp.vcf "
            "--filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0\" "
            "--filter-name \"basic_snp_filter\" "
            "-O {params.pre}/03.snp.filtered.vcf"),
        shell("{params.cmd} VariantFiltration -R {params.ref} "
            "-V {params.pre}/02.indel.vcf "
            "--filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\" "
            "--filter-name \"basic_indel_filter\" "
            "-O {params.pre}/03.indel.filtered.vcf"),
        shell("{params.cmd} MergeVcfs "
            "-I {params.pre}/03.snp.filtered.vcf "
            "-I {params.pre}/03.indel.filtered.vcf "
            "-O {output}"),
        shell("vcf 2tsv {params.pre}/10.vcf > {params.pre}/10.tsv")


