def gatk_hc_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['haplotype_caller']:
        java_mem = config['gatk']['haplotype_caller']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
    java_options = "-Xmx%s -Djava.io.tmpdir=%s" % (java_mem, java_tmpdir)
    cmd = "%s --java-options \"%s\"" % (config['gatk']['cmd'], java_options)
    return cmd

rule gatk_haplotype_caller:
    input:
        "%s/{sid}.bam" % config['callvnt']['idir']
    output:
        temp("%s/{sid}/{rid}.g.vcf.gz" % config['callvnt']['odir']),
        temp("%s/{sid}/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir'])
    log:
        "%s/gatk/{sid}/{rid}.log" % config['dirl']
    params:
        cmd = gatk_hc_cmd,
        ref = config['gatk']['ref'],
        region = lambda wildcards: config['regions'][wildcards.rid],
        extra = config['gatk']['haplotype_caller']['extra'],
        ppn = config['gatk']['haplotype_caller']['ppn'],
        walltime = config['gatk']['haplotype_caller']['walltime'],
        mem = config['gatk']['haplotype_caller']['mem']
    threads: config['gatk']['haplotype_caller']['ppn']
    shell:
        #-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        """
        {params.cmd} HaplotypeCaller \
        -R {params.ref} \
        -ERC GVCF \
        -L {params.region} \
        {params.extra} \
        -I {input} \
        -O {output[0]} \
        >{log} 2>&1
        """

def gather_vcf_inputs(wildcards):
    sid = wildcards.sid
    vcfs = expand(["%s/%s/{rid}.g.vcf.gz" % (config['callvnt']['odir1'], sid)],
                  rid = list(config['regions'].keys()))
    tbis = expand(["%s/%s/{rid}.g.vcf.gz.tbi" % (config['callvnt']['odir1'], sid)],
                  rid = list(config['regions'].keys()))
    return {
        'vcfs': vcfs,
        'tbis': tbis
    }

rule gatk_gather_vcfs:
    input:
        unpack(gather_vcf_inputs)
    output:
        vcf = protected("%s/{sid}.g.vcf.gz" % config['callvnt']['odir1']),
        tbi = protected("%s/{sid}.g.vcf.gz.tbi" % config['callvnt']['odir1'])
    log:
        "%s/gatk/{sid}.log" % config['dirl']
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda wildcards, input: ["-I %s" % x for x in input.vcfs],
        extra = '',
        ppn = config['gatk']['gather_vcfs']['ppn'],
        walltime = config['gatk']['gather_vcfs']['walltime'],
        mem = config['gatk']['gather_vcfs']['mem']
    threads: config['gatk']['gather_vcfs']['ppn']
    shell:
        """
        {params.cmd} GatherVcfs {params.input_str} -O {output.vcf} >{log} 2>&1
        tabix -p vcf {output.vcf}
        """

rule gatk_combine_gvcfs:
    input:
        vcfs = expand(["%s/{sid}.g.vcf.gz" % config['callvnt']['odir1']], 
                sid = config['SampleID']),
        tbis = expand(["%s/{sid}.g.vcf.gz.tbi" % config['callvnt']['odir1']], 
                sid = config['SampleID'])
    output:
        vcf = protected("%s/all.g.vcf.gz" % config['callvnt']['odir2']),
        tbi = protected("%s/all.g.vcf.gz.tbi" % config['callvnt']['odir2']),
    params:
        cmd = gatk_gg_cmd,
        ref = config['gatk']['ref'],
        gvcfs = lambda wildcards, input: ["-V %s" % x for x in input.vcfs],
        ppn = config['gatk']['combine_gvcfs']['ppn'],
    threads: config['gatk']['combine_gvcfs']['ppn']
    shell:
        """
        {params.cmd} CombineGVCFs \
        -R {params.ref} \
        {params.gvcfs} -O {output.vcf}
        """

def gatk_gg_cmd(wildcards):
    java_mem = config['java']['mem']
    if 'java_mem' in config['gatk']:
        java_mem = config['gatk']['java_mem']
    if 'java_mem' in config['gatk']['genotype_gvcfs']:
        java_mem = config['gatk']['genotype_gvcfs']['java_mem']
    java_tmpdir = config['tmpdir']
    if 'tmpdir' in config['java']:
        java_tmpdir = config['java']['tmpdir']
    java_options = "-Xmx%s -Djava.io.tmpdir=%s" % (java_mem, java_tmpdir)
    cmd = "%s --java-options \"%s\"" % (config['gatk']['cmd'], java_options)
    return cmd

rule gatk_genotype_gvcfs:
    input:
        vcf = "%s/all.g.vcf.gz" % config['callvnt']['odir2'],
        tbi = "%s/all.g.vcf.gz.tbi" % config['callvnt']['odir2'],
    output:
        vcf = protected("%s/{rid}.vcf.gz" % config['callvnt']['odir2']),
        tbi = protected("%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2']),
    params:
        cmd = gatk_gg_cmd,
        ref = config['gatk']['ref'],
        region = lambda wildcards: config['regions'][wildcards.rid],
        ppn = config['gatk']['genotype_gvcfs']['ppn'],
        walltime = config['gatk']['genotype_gvcfs']['walltime'],
        mem = config['gatk']['genotype_gvcfs']['mem']
    threads: config['gatk']['genotype_gvcfs']['ppn']
    shell:
        """
        {params.cmd} GenotypeGVCFs \
        -nt {threads} \
        -R {params.ref} -L {params.region} \
        -V {input.vcf} -O {output.vcf}
        """

rule gatk_gather_vcfs2:
    input:
        vcfs = expand(["%s/{rid}.vcf.gz" % config['callvnt']['odir2']], 
                rid = list(config['regions'].keys())),
        tbis = expand(["%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2']], 
                rid = list(config['regions'].keys()))
    output:
        vcf = protected("%s" % config['callvnt']['outfile']),
        tbi = protected("%s.tbi" % config['callvnt']['outfile'])
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda wildcards, input: ["-I %s" % x for x in input.vcfs],
    threads:
        config["gatk"]['gather_vcfs']["threads"]
    shell:
        """
        {params.cmd} GatherVcfs {params.input_str} -O {output.vcf}
        tabix -p vcf {output.vcf}
        """
