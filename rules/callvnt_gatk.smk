def gatk_hc_inputs(w):
    sids = config['gt'][w.gt]
    return expand("%s/{sid}.bam" % config['callvnt']['idir'], sid = sids)

rule gatk_haplotype_caller:
    input:
        gatk_hc_inputs
    output:
        temp("%s/{gt}/{rid}.g.vcf.gz" % config['callvnt']['odir1']),
        temp("%s/{gt}/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir1'])
    log:
        "%s/%s/{gt}/{rid}.log" % (config['dirl'], config['gatk']['haplotype_caller']['id'])
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        input_str = lambda w, input: ["-I %s" % x for x in input],
        region = lambda w: config[config['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True, hc = True),
        N = lambda w: "%s.%s.%s" % (config['gatk']['haplotype_caller']['id'], w.gt, w.rid),
        e = lambda w: "%s/%s/%s/%s.e" % (config['dirp'], config['gatk']['haplotype_caller']['id'], w.gt, w.rid),
        o = lambda w: "%s/%s/%s/%s.o" % (config['dirp'], config['gatk']['haplotype_caller']['id'], w.gt, w.rid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem - 2
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['mem'],
        load = lambda w, attempt:  get_resource(config, attempt, 'gatk','haplotype_caller')['load']
    threads: config['gatk']['haplotype_caller']['ppn']
    shell:
        #-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" HaplotypeCaller \
        -R {params.ref} \
        -ERC GVCF \
        -L {params.region} \
        {params.extra} \
        {params.input_str} -O {output[0]} \
        >>{log} 2>&1
        """

def gather_vcf_inputs(w):
    gt = w.gt
    vcfs = expand("%s/%s/{rid}.g.vcf.gz" % (config['callvnt']['odir1'], gt),
                  rid = list(config[config['reference']]['regions'].keys()))
    tbis = expand("%s/%s/{rid}.g.vcf.gz.tbi" % (config['callvnt']['odir1'], gt),
                  rid = list(config[config['reference']]['regions'].keys()))
    return {
        'vcfs': vcfs,
        'tbis': tbis
    }

rule gatk_gather_vcfs:
    input:
        unpack(gather_vcf_inputs)
    output:
        vcf = protected("%s/%s/{gt}.g.vcf.gz" % (config['dird'], config['callvnt']['odir1'])),
        tbi = protected("%s/%s/{gt}.g.vcf.gz.tbi" % (config['dird'], config['callvnt']['odir1']))
    log:
        "%s/%s/{gt}.log" % (config['dirl'], config['gatk']['gather_vcfs']['id'])
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True),
        N = lambda w: "%s.%s" % (config['gatk']['gather_vcfs']['id'], w.gt),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['gather_vcfs']['id'], w.gt),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['gather_vcfs']['id'], w.gt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'gather_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'gather_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'gather_vcfs')['mem']
    threads: config['gatk']['gather_vcfs']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GatherVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf} >{log} 2>&1
        tabix -p vcf {output.vcf}
        """

rule gatk_combine_gvcfs:
    input:
        vcfs = expand("%s/%s/{gt}.g.vcf.gz" %
            (config['dird'], config['callvnt']['odir1']),
            gt = config['Genotypes']),
        tbis = expand("%s/%s/{gt}.g.vcf.gz.tbi" %
            (config['dird'], config['callvnt']['odir1']),
            gt = config['Genotypes'])
    output:
        vcf = "%s/{rid}.g.vcf.gz" % config['callvnt']['odir2'],
        tbi = "%s/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir2'],
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        region = lambda w: config[config['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['combine_gvcfs']['id'], w.rid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['combine_gvcfs']['id'], w.rid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['combine_gvcfs']['id'], w.rid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['mem']
    threads: config['gatk']['combine_gvcfs']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" CombineGVCFs \
        -R {params.ref} -L {params.region} \
        {params.extra} \
        {params.gvcfs} -O {output.vcf}
        """

rule gatk_genotype_gvcfs:
    input:
        vcf = "%s/{rid}.g.vcf.gz" % config['callvnt']['odir2'],
        tbi = "%s/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir2'],
    output:
        vcf = "%s/{rid}.vcf.gz" % config['callvnt']['odir2'],
        tbi = "%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2'],
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        region = lambda w: config[config['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['genotype_gvcfs']['id'], w.rid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['genotype_gvcfs']['id'], w.rid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['genotype_gvcfs']['id'], w.rid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['mem']
    threads: config['gatk']['genotype_gvcfs']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GenotypeGVCFs \
        -R {params.ref} \
        {params.extra} \
        -V {input.vcf} -O {output.vcf}
        """

rule gatk_gather_vcfs2:
    input:
        vcfs = expand("%s/{rid}.vcf.gz" % config['callvnt']['odir2'], 
                rid = list(config[config['reference']]['regions'].keys())),
        tbis = expand("%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2'], 
                rid = list(config[config['reference']]['regions'].keys()))
    output:
        vcf = protected("%s/%s" % (config['dird'], config['callvnt']['out'])),
        tbi = protected("%s/%s.tbi" % (config['dird'], config['callvnt']['out']))
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True),
        N = lambda w: "%s" % (config['gatk']['gather_vcfs']['id']),
        e = lambda w: "%s/%s.e" % (config['dirp'], config['gatk']['gather_vcfs']['id']),
        o = lambda w: "%s/%s.o" % (config['dirp'], config['gatk']['gather_vcfs']['id']),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'gather_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'gather_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'gather_vcfs')['mem']
    threads: config["gatk"]['gather_vcfs']["ppn"]
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GatherVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf}
        bcftools index -t {output.vcf}
        """
