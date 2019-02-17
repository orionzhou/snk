def gatk_hc_inputs(w):
    yid, gt = w.yid, w.gt
    sids = config['y'][yid]['gt'][gt]
    return expand("%s/%s/{sid}.bam" % (yid, config['cleanbam']['odir2']), sid = sids)

rule gatk_haplotype_caller:
    input:
        gatk_hc_inputs
    output:
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz" % config['callvnt']['odir']),
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir'])
    log:
        "{yid}/%s/%s/{gt}/{rid}.log" % (config['dirl'], config['gatk']['haplotype_caller']['id'])
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        input_str = lambda w, input: ["-I %s" % x for x in input],
        region = lambda w: config[config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True, hc = True),
        N = "{yid}.%s.{gt}.{rid}" % config['gatk']['haplotype_caller']['id'],
        e = "{yid}/%s/%s/{gt}/{rid}.e" % (config['dirp'], config['gatk']['haplotype_caller']['id']),
        o = "{yid}/%s/%s/{gt}/{rid}.o" % (config['dirp'], config['gatk']['haplotype_caller']['id']),
        mem = lambda w, resources: resources.mem
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
        {params.extra} -ERC GVCF \
        -R {params.ref} \
        -L {params.region} \
        {params.input_str} -O {output[0]} \
        >>{log} 2>&1
        """

def merge_vcf_inputs(w):
    yid, gt = w.yid, w.gt
    rids= list(config[config['y'][yid]['reference']]['regions'].keys())
    vcfs = expand("%s/%s/%s/{rid}.g.vcf.gz" %
        (yid, config['callvnt']['odir'], gt), rid = rids)
    tbis = expand("%s/%s/%s/{rid}.g.vcf.gz.tbi" %
        (yid, config['callvnt']['odir'], gt), rid = rids)
    return {'vcfs':vcfs, 'tbis':tbis}

rule gatk_merge_vcfs:
    input:
        unpack(merge_vcf_inputs)
    output:
        vcf = protected("{yid}/%s/%s/{gt}.g.vcf.gz" % (config['dird'], config['callvnt']['odir'])),
        tbi = protected("{yid}/%s/%s/{gt}.g.vcf.gz.tbi" % (config['dird'], config['callvnt']['odir']))
    log:
        "{yid}/%s/%s/{gt}.log" % (config['dirl'], config['gatk']['merge_vcfs']['id'])
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True),
        N = "{yid}.%s.{gt}" % config['gatk']['merge_vcfs']['id'],
        e = "{yid}/%s/%s/{gt}.e" % (config['dirp'], config['gatk']['merge_vcfs']['id']),
        o = "{yid}/%s/%s/{gt}.o" % (config['dirp'], config['gatk']['merge_vcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['mem']
    threads: config['gatk']['merge_vcfs']['ppn']
    shell:
        #tabix -p vcf {output.vcf}
        #bcftools index -t {output.vcf}
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf} >{log} 2>&1
        """

#-----------
def gatk_rename_sample_id_inputs(w):
    yid, sid = w.yid, w.sid
    study = config['y'][yid]['t']['yid']
    gt = config['y'][yid]['t']['Genotype']
    vcf = "%s/data/%s/%s.g.vcf.gz" % (study, config['callvnt']['odir'], gt)
    tbi = "%s/data/%s/%s.g.vcf.gz.tbi" % (study, config['callvnt']['odir'], gt)
    return {'vcf':vcf, 'tbi':tbi}

rule gatk_rename_sample_id:
    input: gatk_rename_sample_id_inputs
    output:
        vcf = "{yid}/%s/{sid}.g.vcf.gz" % (config['callvnt']['odir1']),
        tbi = "{yid}/%s/{sid}.g.vcf.gz.tbi" % (config['callvnt']['odir1']),
    params:
        cmd = config['gatk']['cmd'],
        extra = gatk_extra(picard = True, jdk = False),
        ngt = lambda w: "%s_%s" % (w.study, w.gt),
        N = lambda w: "%s.%s.%s" % (config['gatk']['rename_sample_id']['id'], w.study,w.gt),
        e = lambda w: "%s/%s/%s/%s.e" % (config['dirp'], config['gatk']['rename_sample_id']['id'], w.study,w.gt),
        o = lambda w: "%s/%s/%s/%s.o" % (config['dirp'], config['gatk']['rename_sample_id']['id'], w.study,w.gt),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['mem']
    threads: config['gatk']['rename_sample_id']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" RenameSampleInVcf \
        {params.extra} \
        --NEW_SAMPLE_NAME {params.ngt} --OLD_SAMPLE_NAME {wildcards.gt} \
        -I {input.vcf} -O {output.vcf} --CREATE_INDEX
        """

def combine_gvcf_inputs(w):
    yid = w.yid
    vcfs, tbis = [], []
    for sid, sdic in config['y'][yid]['t'].items():
        fv = '%s/%s_%s.g.vcf.gz' % (config['callvnt']['odir1'], sid)
        fx = '%s/%s_%s.g.vcf.gz.tbi' % (config['callvnt']['odir1'], sid)
        vcfs.append(fv)
        tbis.append(fx)
    return {'vcfs':vcfs, 'tbis':tbis}

rule gatk_combine_gvcfs:
    input:
        vcfs = lambda w: expand("%s/%s/{sid}.g.vcf.gz" % (w.yid, config['callvnt']['odir1']), sid = config['y'][w.yid]['t'].keys()),
        tbis = lambda w: expand("%s/%s/{sid}.g.vcf.gz.tbi" % (w.yid, config['callvnt']['odir1']), sid = config['y'][w.yid]['t'].keys()),
    output:
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt']['odir2'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir2'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        region = lambda w: config[config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['combine_gvcfs']['id'], w.rid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['combine_gvcfs']['id'], w.rid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['combine_gvcfs']['id'], w.rid),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['mem']
    threads: config['gatk']['combine_gvcfs']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" CombineGVCFs \
        {params.extra} \
        -R {params.ref} -L {params.region} \
        {params.gvcfs} -O {output.vcf}
        """

rule gatk_genotype_gvcfs:
    input:
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt']['odir2'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt']['odir2'],
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt']['odir3'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir3'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        region = lambda w: config[config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['genotype_gvcfs']['id'], w.rid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['genotype_gvcfs']['id'], w.rid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['genotype_gvcfs']['id'], w.rid),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['mem']
    threads: config['gatk']['genotype_gvcfs']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GenotypeGVCFs \
        {params.extra} \
        -R {params.ref} \
        -V {input.vcf} -O {output.vcf}
        """

rule gatk_merge_vcfs2:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.vcf.gz" % (w.yid, config['callvnt']['odir3']),
                rid = config[config['y'][w.yid]['reference']]['regions'].keys()),
        tbis = lambda w: expand("%s/%s/{rid}.vcf.gz.tbi" % (w.yid, config['callvnt']['odir3']),
                rid = config[config['y'][w.yid]['reference']]['regions'].keys()),
    output:
        vcf = protected("{yid}/%s/%s" % (config['dird'], config['callvnt']['out'])),
        tbi = protected("{yid}/%s/%s.tbi" % (config['dird'], config['callvnt']['out']))
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True, jdk = False),
        N = lambda w: "%s" % (config['gatk']['merge_vcfs']['id']),
        e = lambda w: "%s/%s.e" % (config['dirp'], config['gatk']['merge_vcfs']['id']),
        o = lambda w: "%s/%s.o" % (config['dirp'], config['gatk']['merge_vcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['mem']
    threads: config["gatk"]['merge_vcfs']["ppn"]
    shell:
        #bcftools index -t {output.vcf}
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf}
        """
