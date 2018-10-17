def gatk_hc_walltime(attempt):
    seconds = config['gatk']['haplotype_caller']['walltime']
    seconds_extra = (attempt-1) * 15 * 3600 
    return seconds + seconds_extra 

def gatk_hc_mem(attempt):
    memstr = config['gatk']['haplotype_caller']['mem']
    mem_gb = int(memstr.replace("GB", ""))
    return (mem_gb + (attempt-1) * 10)

rule gatk_haplotype_caller:
    input:
        "%s/{sid}.bam" % config['callvnt']['idir']
    output:
        temp("%s/{sid}/{rid}.g.vcf.gz" % config['callvnt']['odir1']),
        temp("%s/{sid}/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir1'])
    log:
        "%s/gatk/{sid}/{rid}.log" % config['dirl']
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        region = lambda w: config[config['reference']]['regions'][w.rid],
        tmp = config['tmpdir'],
        extra = '-pairHMM LOGLESS_CACHING --use-jdk-deflater --use-jdk-inflater',
        N = lambda w: "hc.%s.%s" % (w.sid, w.rid),
        ppn = config['gatk']['haplotype_caller']['ppn'],
        #walltime = config['gatk']['haplotype_caller']['walltime'],
        #mem = config['gatk']['haplotype_caller']['mem'],
        walltime = lambda w, resources: resources.walltime,
        mem = lambda w, resources: '%dGB' % resources.mem
    resources:
        walltime = lambda w, attempt: gatk_hc_walltime(attempt),
        mem = lambda w, attempt: gatk_hc_mem(attempt),
    threads: config['gatk']['haplotype_caller']['ppn']
    shell:
        #-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        """
        {params.cmd} \
        --java_options "-Xmx{params.mem} -Djava.io.tmpdir={params.tmp}" \
        HaplotypeCaller \
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
                  rid = list(config[config['reference']]['regions'].keys()))
    tbis = expand(["%s/%s/{rid}.g.vcf.gz.tbi" % (config['callvnt']['odir1'], sid)],
                  rid = list(config[config['reference']]['regions'].keys()))
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
        tmp = config['tmpdir'],
        N = lambda w: "gather.%s" % (w.sid),
        ppn = config['gatk']['gather_vcfs']['ppn'],
        walltime = config['gatk']['gather_vcfs']['walltime'],
        mem = config['gatk']['gather_vcfs']['mem']
    threads: config['gatk']['gather_vcfs']['ppn']
    shell:
        """
        {params.cmd} \
        --java_options "-Xmx{params.mem} -Djava.io.tmpdir={params.tmp}" \
        GatherVcfs \
        {params.input_str} -O {output.vcf} >{log} 2>&1
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
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        gvcfs = lambda wildcards, input: ["-V %s" % x for x in input.vcfs],
        extra = '-pairHMM LOGLESS_CACHING --use-jdk-deflater --use-jdk-inflater',
        tmp = config['tmpdir'],
        N = lambda w: 'combine.gvcf',
        ppn = config['gatk']['combine_gvcfs']['ppn'],
        walltime = config['gatk']['combine_gvcfs']['walltime'],
        mem = config['gatk']['combine_gvcfs']['mem']
    threads: config['gatk']['combine_gvcfs']['ppn']
    shell:
        """
        {params.cmd} \
        --java_options "-Xmx{params.mem} -Djava.io.tmpdir={params.tmp}" \
        CombineGVCFs \
        -R {params.ref} \
        {params.extra} \
        {params.gvcfs} -O {output.vcf}
        """

rule gatk_genotype_gvcfs:
    input:
        vcf = "%s/all.g.vcf.gz" % config['callvnt']['odir2'],
        tbi = "%s/all.g.vcf.gz.tbi" % config['callvnt']['odir2'],
    output:
        vcf = protected("%s/{rid}.vcf.gz" % config['callvnt']['odir2']),
        tbi = protected("%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2']),
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        region = lambda wildcards: config[config['reference']]['regions'][wildcards.rid],
        extra = '-pairHMM LOGLESS_CACHING --use-jdk-deflater --use-jdk-inflater',
        tmp = config['tmpdir'],
        N = lambda w: 'gg.%s' % w.rid,
        ppn = config['gatk']['genotype_gvcfs']['ppn'],
        walltime = config['gatk']['genotype_gvcfs']['walltime'],
        mem = config['gatk']['genotype_gvcfs']['mem']
    threads: config['gatk']['genotype_gvcfs']['ppn']
    shell:
        """
        {params.cmd} \
        --java_options "-Xmx{params.mem} -Djava.io.tmpdir={params.tmp}" \
        GenotypeGVCFs \
        -nt {threads} \
        {params.extra} \
        -R {params.ref} -L {params.region} \
        -V {input.vcf} -O {output.vcf}
        """

rule gatk_gather_vcfs2:
    input:
        vcfs = expand(["%s/{rid}.vcf.gz" % config['callvnt']['odir2']], 
                rid = list(config[config['reference']]['regions'].keys())),
        tbis = expand(["%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2']], 
                rid = list(config[config['reference']]['regions'].keys()))
    output:
        vcf = protected("%s" % config['callvnt']['outfile']),
        tbi = protected("%s.tbi" % config['callvnt']['outfile'])
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda wildcards, input: ["-I %s" % x for x in input.vcfs],
        tmp = config['tmpdir'],
        N = lambda w: 'gather.vcf',
        ppn = config['gatk']['gather_vcfs']['ppn'],
        walltime = config['gatk']['gather_vcfs']['walltime'],
        mem = config['gatk']['gather_vcfs']['mem']
    threads: config["gatk"]['gather_vcfs']["ppn"]
    shell:
        """
        {params.cmd} \
        --java_options "-Xmx{params.mem} -Djava.io.tmpdir={params.tmp}" \
        GatherVcfs \
        {params.input_str} -O {output.vcf}
        bcftools index -t {output.vcf}
        """
