def gatk_hc_inputs(w):
    yid, gt = w.yid, w.gt
    sids = config['y'][yid]['gt'][gt]
    return expand("%s/%s/{sid}.bam" % (yid, config['cleanbam']['od24']), sid = sids)

rule gatk_haplotype_caller:
    input: gatk_hc_inputs
    output:
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz" % config['callvnt']['od25']),
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz.tbi" % config['callvnt']['od25'])
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        input_str = lambda w, input: ["-I %s" % x for x in input],
        region = lambda w: config['g'][config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True, hc = True),
        N = "{yid}.%s.{gt}.{rid}" % config['gatk']['haplotype_caller']['id'],
        e = "{yid}/%s/%s/{gt}/{rid}.e" % (config['dirj'], config['gatk']['haplotype_caller']['id']),
        o = "{yid}/%s/%s/{gt}/{rid}.o" % (config['dirj'], config['gatk']['haplotype_caller']['id']),
        mem = lambda w, resources: resources.mem - 3
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'gatk', 'haplotype_caller')['mem'],
        load = lambda w, attempt:  get_resource(config, attempt, 'gatk','haplotype_caller')['load']
    threads: config['gatk']['haplotype_caller']['ppn']
    conda: "../envs/work.yml"
    shell:
        #-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" HaplotypeCaller \
        {params.extra} -ERC GVCF \
        -R {params.ref} \
        -L {params.region} \
        {params.input_str} -O {output[0]}
        """

def merge_vcf_inputs(w):
    yid, gt = w.yid, w.gt
    rids= list(config['g'][config['y'][yid]['reference']]['regions'].keys())
    vcfs = expand("%s/%s/%s/{rid}.g.vcf.gz" %
        (yid, config['callvnt']['od25'], gt), rid = rids)
    tbis = expand("%s/%s/%s/{rid}.g.vcf.gz.tbi" %
        (yid, config['callvnt']['od25'], gt), rid = rids)
    return {'vcfs':vcfs, 'tbis':tbis}

rule gatk_merge_vcfs:
    input: unpack(merge_vcf_inputs)
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
        e = "{yid}/%s/%s/{gt}.e" % (config['dirj'], config['gatk']['merge_vcfs']['id']),
        o = "{yid}/%s/%s/{gt}.o" % (config['dirj'], config['gatk']['merge_vcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['mem']
    threads: config['gatk']['merge_vcfs']['ppn']
    conda: "../envs/work.yml"
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
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['gatk']['rename_sample_id']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['gatk']['rename_sample_id']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'rename_sample_id')['mem']
    threads: config['gatk']['rename_sample_id']['ppn']
    conda: "../envs/work.yml"
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
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        region = lambda w: config['g'][config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk']['combine_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['gatk']['combine_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['gatk']['combine_gvcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'combine_gvcfs')['mem']
    threads: config['gatk']['combine_gvcfs']['ppn']
    conda: "../envs/work.yml"
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
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        region = lambda w: config['g'][config['y'][w.yid]['reference']]['regions'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk']['genotype_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['gatk']['genotype_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['gatk']['genotype_gvcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'genotype_gvcfs')['mem']
    threads: config['gatk']['genotype_gvcfs']['ppn']
    conda: "../envs/work.yml"
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
                rid = config['g'][config['y'][w.yid]['reference']]['regions'].keys()),
        tbis = lambda w: expand("%s/%s/{rid}.vcf.gz.tbi" % (w.yid, config['callvnt']['od28']),
                rid = config['g'][config['y'][w.yid]['reference']]['regions'].keys()),
    output:
        vcf = "{yid}/%s" % config['callvnt']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of30']
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True, jdk = False),
        N = "{yid}.%s" % config['gatk']['merge_vcfs']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk']['merge_vcfs']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk']['merge_vcfs']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'merge_vcfs')['mem']
    threads: config["gatk"]['merge_vcfs']["ppn"]
    conda: "../envs/work.yml"
    shell:
        #bcftools index -t {output.vcf}
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf}
        """

#___ variant recal
rule gatk_pick_training:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of30'],
    output:
        snp = "{yid}/%s" % config['callvnt']['of31a'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt']['of31a'],
        idl = "{yid}/%s" % config['callvnt']['of31b'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt']['of31b'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        minqual = 990,
        N = "%s.{yid}" % config['gatk']['pick_training']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk']['pick_training']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk']['pick_training']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'pick_training')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'pick_training')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'pick_training')['mem']
    threads: config["gatk"]['pick_training']["ppn"]
    conda: "../envs/work.yml"
    shell:
        """
        bcftools view -G -i 'TYPE="snp" & QUAL>{params.minqual}' {input.vcf} -Oz -o {output.snp}
        bcftools view -G -i 'TYPE="indel" & QUAL>{params.minqual}' {input.vcf} -Oz -o {output.idl}
        bcftools index -t {output.snp}
        bcftools index -t {output.idl}
        bcftools stats {output.snp} > {output.snp}.txt
        bcftools stats {output.idl} > {output.idl}.txt
        """

rule gatk_variant_recalibrator:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of30'],
#        truth = config['callvnt']['truth_sites_snp'],
#        truth_tbi = config['callvnt']['truth_sites_snp'].replace('.gz','.gz.tbi'),
        snp = "{yid}/%s" % config['callvnt']['of31a'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt']['of31a'],
        idl = "{yid}/%s" % config['callvnt']['of31b'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt']['of31b'],
    output:
        snp_vcf = "{yid}/%s" % config['callvnt']['of33'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt']['of33'],
        idl_vcf = "{yid}/%s" % config['callvnt']['of35'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt']['of35']
    params:
        snp_recal = "{yid}/32.recal.vcf",
        snp_tranch = "{yid}/32.tranches",
        snp_model = "{yid}/32.model",
        idl_recal = "{yid}/34.recal.vcf",
        idl_tranch = "{yid}/34.tranches",
        idl_model = "{yid}/34.model",
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "%s.{yid}" % config['gatk']['variant_recalibrator']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk']['variant_recalibrator']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk']['variant_recalibrator']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_recalibrator')['mem']
    threads: config["gatk"]['variant_recalibrator']["ppn"]
    conda: "../envs/work.yml"
    shell:
# 	--output-model {output.report} \
#       --maxGaussians 4 \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" VariantRecalibrator \
        {params.extra} -R {params.ref} \
	--resource:hmp3,known=false,training=true,truth=true,prior=10.0 {input.snp} \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
 	-mode SNP \
        -V {input.vcf} -O {params.snp_recal} --tranches-file {params.snp_tranch}

        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyVQSR \
        {params.extra} -R {params.ref} \
        -V {input.vcf} -O {output.snp_vcf} \
	--recal-file {params.snp_recal} --tranches-file {params.snp_tranch} \
	--truth-sensitivity-filter-level 99.5 \
 	-mode SNP \
 	--create-output-variant-index true

        {params.cmd} --java-options "-Xmx{params.mem}G" VariantRecalibrator \
        {params.extra} -R {params.ref} \
	--resource:hmp3,known=false,training=true,truth=true,prior=10.0 {input.idl} \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
        -mode INDEL \
        -V {output.snp_vcf} -O {params.idl_recal} --tranches-file {params.idl_tranch}

        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyVQSR \
        {params.extra} -R {params.ref} \
        -V {output.snp_vcf} -O {output.idl_vcf} \
	--recal-file {params.idl_recal} --tranches-file {params.idl_tranch} \
	--truth-sensitivity-filter-level 99.0 \
	-mode INDEL \
 	--create-output-variant-index true
        """

# final filter
rule gatk_variant_filtration:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of35'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of35'],
    output:
        vcf = "{yid}/%s" % config['callvnt']['of37'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of37'],
        stat = "{yid}/%s" % config['callvnt']['of37'].replace(".vcf.gz", ".txt")
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        N = "%s.{yid}" % config['gatk']['variant_filtration']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk']['variant_filtration']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk']['variant_filtration']['id']),
        mem = lambda w, resources: resources.mem - 1
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_filtration')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_filtration')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'variant_filtration')['mem']
    threads: config["gatk"]['variant_filtration']["ppn"]
    conda: "../envs/work.yml"
    shell:
        """
        bcftools filter -i 'FILTER=="PASS"' {input.vcf} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        bcftools stats -s - {output.vcf} > {output.stat}
        """

rule gatk_snpeff:
    input:
        vcf = "{yid}/%s" % config['callvnt']['of37'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of37'],
    output:
        vcf = "{yid}/%s" % config['callvnt']['of38'],
        tbi = "{yid}/%s.tbi" % config['callvnt']['of38'],
    params:
        genome = lambda w: config['y'][w.yid]['reference'],
        idx = lambda w: config['g'][config['y'][w.yid]['reference']]['snpeff']['xcfg'],
        N = "%s.{yid}" % config['snpeff']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['snpeff']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['snpeff']['id']),
        mem = lambda w, resources: resources.mem - 1
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'snpeff')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'snpeff')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'snpeff')['mem']
    threads: config['snpeff']["ppn"]
    conda: "../envs/work.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.idx} {params.genome} -ud 0 \
                {input.vcf} | bgzip > {output.vcf}
        bcftools index -t {output.vcf}
        """


