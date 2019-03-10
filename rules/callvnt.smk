def gatk_hc_inputs(w):
    yid, gt = w.yid, w.gt
    sids = config['y'][yid]['gt'][gt]
    return expand("%s/%s/{sid}.bam" % (yid, config['cleanbam']['odir2']), sid = sids)

rule gatk_haplotype_caller:
    input: gatk_hc_inputs
    output:
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz" % config['callvnt']['od25']),
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz.tbi" % config['callvnt']['od25'])
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
        mem = lambda w, resources: resources.mem - 3
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['mem'],
        load = lambda w, attempt:  get_resource(config, attempt, 'gatk','haplotype_caller')['load']
    threads: config['gatk']['haplotype_caller']['ppn']
    conda: "../envs/gatk.yml"
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
        (yid, config['callvnt']['od25'], gt), rid = rids)
    tbis = expand("%s/%s/%s/{rid}.g.vcf.gz.tbi" %
        (yid, config['callvnt']['od25'], gt), rid = rids)
    return {'vcfs':vcfs, 'tbis':tbis}

rule gatk_merge_vcfs:
    input:
        unpack(merge_vcf_inputs)
    output:
        vcf = protected("%s/{yid}/{gt}.g.vcf.gz" % config['callvnt']['odir']),
        tbi = protected("%s/{yid}/{gt}.g.vcf.gz.tbi" % config['callvnt']['odir'])
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
    conda: "../envs/gatk.yml"
    shell:
        #tabix -p vcf {output.vcf}
        #bcftools index -t {output.vcf}
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf} >{log} 2>&1
        """

#----
def gatk_rename_sample_id_inputs(w):
    yid, sid = w.yid, w.sid
    study = config['y'][yid]['t'][sid]['yid']
    gt = config['y'][yid]['t'][sid]['Genotype']
    vcf = "%s/%s/%s.g.vcf.gz" % (config['callvnt']['odir'], study, gt)
    tbi = "%s/%s/%s.g.vcf.gz.tbi" % (config['callvnt']['odir'], study, gt)
    return {'vcf':vcf, 'tbi':tbi}

rule gatk_rename_sample_id:
    input: unpack(gatk_rename_sample_id_inputs)
    output:
        vcf = "{yid}/%s/{sid}.g.vcf.gz" % config['callvnt']['od26'],
        tbi = "{yid}/%s/{sid}.g.vcf.gz.tbi" % config['callvnt']['od26'],
    params:
        cmd = config['gatk']['cmd'],
        extra = gatk_extra(picard = True, jdk = False),
        ogt = lambda w: config['y'][w.yid]['t'][w.sid]['Genotype'],
        N = "{yid}.%s.{sid}" % config['gatk']['rename_sample_id']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['gatk']['rename_sample_id']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['gatk']['rename_sample_id']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['mem']
    threads: config['gatk']['rename_sample_id']['ppn']
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" RenameSampleInVcf \
        {params.extra} \
        --NEW_SAMPLE_NAME {wildcards.sid} --OLD_SAMPLE_NAME {params.ogt} \
        -I {input.vcf} -O {output.vcf} --CREATE_INDEX
        """

def combine_gvcf_inputs(w):
    yid = w.yid
    vcfs, tbis = [], []
    for sid, sdic in config['y'][yid]['t'].items():
        fv = '%s/%s.g.vcf.gz' % (config['callvnt']['od26'], sid)
        fx = '%s/%s.g.vcf.gz.tbi' % (config['callvnt']['od26'], sid)
        vcfs.append(fv)
        tbis.append(fx)
    return {'vcfs':vcfs, 'tbis':tbis}

rule gatk_combine_gvcfs:
    input:
        vcfs = lambda w: expand("%s/%s/{sid}.g.vcf.gz" % (w.yid, config['callvnt']['od26']), sid = config['y'][w.yid]['t'].keys()),
        tbis = lambda w: expand("%s/%s/{sid}.g.vcf.gz.tbi" % (w.yid, config['callvnt']['od26']), sid = config['y'][w.yid]['t'].keys()),
    output:
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt']['od27'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt']['od27'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        region = lambda w: config[config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk']['combine_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirp'], config['gatk']['combine_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirp'], config['gatk']['combine_gvcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['mem']
    threads: config['gatk']['combine_gvcfs']['ppn']
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" CombineGVCFs \
        {params.extra} \
        -R {params.ref} -L {params.region} \
        {params.gvcfs} -O {output.vcf}
        """

rule gatk_genotype_gvcfs:
    input:
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt']['od27'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt']['od27'],
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt']['od28'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt']['od28'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        region = lambda w: config[config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk']['genotype_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirp'], config['gatk']['genotype_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirp'], config['gatk']['genotype_gvcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['mem']
    threads: config['gatk']['genotype_gvcfs']['ppn']
    conda: "../envs/gatk.yml"
    shell:
#        -L {input.bed} --include-non-variant-sites \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GenotypeGVCFs \
        {params.extra} \
        -R {params.ref} \
        -V {input.vcf} -O {output.vcf}
        """

rule gatk_merge_vcfs2:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.vcf.gz" % (w.yid, config['callvnt']['od28']),
                rid = config[config['y'][w.yid]['reference']]['regions'].keys()),
        tbis = lambda w: expand("%s/%s/{rid}.vcf.gz.tbi" % (w.yid, config['callvnt']['od28']),
                rid = config[config['y'][w.yid]['reference']]['regions'].keys()),
    output:
        vcf = "{yid}/%s" % config['callvnt']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of30']
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True, jdk = False),
        N = "{yid}.%s" % config['gatk']['merge_vcfs']['id'],
        e = "{yid}/%s/%s.e" % (config['dirp'], config['gatk']['merge_vcfs']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['gatk']['merge_vcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['mem']
    threads: config["gatk"]['merge_vcfs']["ppn"]
    conda: "../envs/gatk.yml"
    shell:
        #bcftools index -t {output.vcf}
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf}
        """

#___ variant recal
rule gatk_variant_recalibrator:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of30'],
	truth = config['callvnt']['truth_sites_snp'],
	truth_tbi = config['callvnt']['truth_sites_snp'].replace('.gz','.gz.tbi'),
    output:
        recal = "{yid}/%s" % config['callvnt']['of32'].replace(".vcf.gz", ".recal.vcf"),
        tranch = "{yid}/%s" % config['callvnt']['of32'].replace(".vcf.gz", ".tranches"),
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "%s.{yid}" % config['gatk']['variant_recalibrator']['id'],
        e = "{yid}/%s/%s.e" % (config['dirp'], config['gatk']['variant_recalibrator']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['gatk']['variant_recalibrator']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['mem']
    threads: config["gatk"]['variant_recalibrator']["ppn"]
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" VariantRecalibrator \
        {params.extra} \
	-R {params.ref} \
	--resource:hmp3,known=false,training=true,truth=true,prior=10.0 {input.truth} \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
 	-mode SNP \
        -V {input.vcf} -O {output.recal} --tranches-file {output.tranch}
        """

rule gatk_apply_vqsr:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of30'],
        recal = "{yid}/%s" % config['callvnt']['of32'].replace(".vcf.gz", ".recal.vcf"),
        tranch = "{yid}/%s" % config['callvnt']['of32'].replace(".vcf.gz", ".tranches"),
    output:
        vcf = "{yid}/%s" % config['callvnt']['of32'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of32']
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(jdk = False),
        N = "%s.{yid}" % config['gatk']['apply_vqsr']['id'],
        e = "{yid}/%s/%s.e" % (config['dirp'], config['gatk']['apply_vqsr']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['gatk']['apply_vqsr']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_vqsr')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_vqsr')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_vqsr')['mem']
    threads: config["gatk"]['apply_vqsr']["ppn"]
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyVQSR \
        {params.extra} \
	-R {params.ref} \
        -V {input.vcf} -O {output.vcf} \
	--tranches-file {input.tranch} --recal-file {input.recal} \
	--truth-sensitivity-filter-level 99.5 \
 	-mode SNP \
 	--create-output-variant-index true
        """

rule gatk_variant_recalibrator2:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of32'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of32'],
	truth = config['callvnt']['truth_sites_idl'],
	truth_tbi = config['callvnt']['truth_sites_idl'].replace('.gz','.gz.tbi'),
    output:
        recal = "{yid}/%s" % config['callvnt']['of34'].replace(".vcf.gz", ".recal.vcf"),
        tranch = "{yid}/%s" % config['callvnt']['of34'].replace(".vcf.gz", ".tranches"),
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(jdk = False),
        N = "%s.{yid}" % config['gatk']['variant_recalibrator']['id'],
        e = "{yid}/%s/%s.e" % (config['dirp'], config['gatk']['variant_recalibrator']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['gatk']['variant_recalibrator']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['mem']
    threads: config["gatk"]['variant_recalibrator']["ppn"]
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" VariantRecalibrator \
        {params.extra} \
	-R {params.ref} \
	--maxGaussians 4 \
	--resource:hmp3,known=false,training=true,truth=true,prior=10.0 {input.truth} \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
        -mode INDEL \
        -V {input.vcf} -O {output.recal} --tranches-file {output.tranch}
        """

rule gatk_apply_vqsr2:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of32'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of32'],
        recal = "{yid}/%s" % config['callvnt']['of34'].replace(".vcf.gz", ".recal.vcf"),
        tranch = "{yid}/%s" % config['callvnt']['of34'].replace(".vcf.gz", ".tranches"),
    output:
        vcf = "{yid}/%s" % config['callvnt']['of34'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of34']
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(jdk = False),
        N = "%s.{yid}" % config['gatk']['apply_vqsr']['id'],
        e = "{yid}/%s/%s.e" % (config['dirp'], config['gatk']['apply_vqsr']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['gatk']['apply_vqsr']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_vqsr')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_vqsr')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_vqsr')['mem']
    threads: config["gatk"]['apply_vqsr']["ppn"]
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyVQSR \
        {params.extra} \
	-R {params.ref} \
        -V {input.vcf} -O {output.vcf} \
	--tranches-file {input.tranch} --recal-file {input.recal} \
	--truth-sensitivity-filter-level 99.0 \
 	-mode INDEL \
 	--create-output-variant-index true
        """




