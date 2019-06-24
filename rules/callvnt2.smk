include: 'gatk.smk'

#def genomics_dbimport_inputs(w):
def combine_gvcfs_inputs(w):
    yid = w.yid
    vcfs, tbis = [], []
    for sid, sdic in config['y'][yid]['t'].items():
        oyid, gt = sdic['yid'], sdic['gt']
        fv = '%s/%s/%s.g.vcf.gz' % (config['callvnt2']['vdir'], oyid, gt)
        fx = '%s/%s/%s.g.vcf.gz.tbi' % (config['callvnt2']['vdir'], oyid, gt)
        vcfs.append(fv)
        tbis.append(fx)
    return {'vcfs':vcfs, 'tbis':tbis}

rule gatk_combine_gvcfs:
    input: unpack(combine_gvcfs_inputs)
    output:
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt2']['od05'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt2']['od05'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win127'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk_combine_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['gatk_combine_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['gatk_combine_gvcfs']['id']),
        j = lambda w: get_resource(w, config, 'gatk_combine_gvcfs'),
        mem = lambda w: get_resource(w, config, 'gatk_combine_gvcfs')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_combine_gvcfs')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" CombineGVCFs \
            {params.extra} -R {params.ref} -L {params.region} \
            {params.gvcfs} -O {output.vcf}
        """

rule gatk_genotype_gvcfs:
    input: #"{yid}/%s/{rid}.cp" % config['callvnt2']['od10']
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt2']['od05'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt2']['od05'],
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt2']['od20'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt2']['od20'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        gendb = "{yid}/%s/{rid}" % config['callvnt2']['od10'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk_genotype_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['gatk_genotype_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['gatk_genotype_gvcfs']['id']),
        j = lambda w: get_resource(w, config, 'gatk_genotype_gvcfs'),
        mem = lambda w: get_resource(w, config, 'gatk_genotype_gvcfs')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_genotype_gvcfs')['ppn']
    conda: "../envs/work.yml"
    shell:
#        -L {input.bed} --include-non-variant-sites \
#        -V gendb://{params.gendb} \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GenotypeGVCFs \
            {params.extra} -R {params.ref} \
            -V {input.vcf} \
            -O {output.vcf}
        """

rule gatk_merge_vcfs2:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.vcf.gz" % (w.yid, config['callvnt2']['od20']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win127'].keys())),
        tbis = lambda w: expand("%s/%s/{rid}.vcf.gz.tbi" % (w.yid, config['callvnt2']['od20']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win127'].keys())),
    output:
        vcf = "{yid}/%s" % config['callvnt2']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of30']
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True, jdk = False),
        N = "{yid}.%s" % config['gatk_merge_vcfs2']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk_merge_vcfs2']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk_merge_vcfs2']['id']),
        j = lambda w: get_resource(w, config, 'gatk_merge_vcfs2'),
        mem = lambda w: get_resource(w, config, 'gatk_merge_vcfs2')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_merge_vcfs2')['ppn']
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
        vcf = "{yid}/%s" % config['callvnt2']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of30'],
    output:
        snp = "{yid}/%s" % config['callvnt2']['of31a'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31a'],
        idl = "{yid}/%s" % config['callvnt2']['of31b'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31b'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        minqual = 990,
        N = "%s.{yid}" % config['gatk_pick_training']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk_pick_training']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk_pick_training']['id']),
        j = lambda w: get_resource(w, config, 'gatk_pick_training'),
        mem = lambda w: get_resource(w, config, 'gatk_pick_training')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_pick_training')['ppn']
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
        vcf = "{yid}/%s" % config['callvnt2']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of30'],
#        truth = config['callvnt']['truth_sites_snp'],
#        truth_tbi = config['callvnt']['truth_sites_snp'].replace('.gz','.gz.tbi'),
        snp = "{yid}/%s" % config['callvnt2']['of31a'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31a'],
        idl = "{yid}/%s" % config['callvnt2']['of31b'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31b'],
    output:
        snp_vcf = "{yid}/%s" % config['callvnt2']['of33'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of33'],
        idl_vcf = "{yid}/%s" % config['callvnt2']['of35'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt2']['of35']
    params:
        snp_recal = "{yid}/32.recal.vcf",
        snp_tranch = "{yid}/32.tranches",
        snp_model = "{yid}/32.model",
        idl_recal = "{yid}/34.recal.vcf",
        idl_tranch = "{yid}/34.tranches",
        idl_model = "{yid}/34.model",
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "%s.{yid}" % config['gatk_variant_recalibrator']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk_variant_recalibrator']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk_variant_recalibrator']['id']),
        j = lambda w: get_resource(w, config, 'gatk_variant_recalibrator'),
        mem = lambda w: get_resource(w, config, 'gatk_variant_recalibrator')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_variant_recalibrator')['ppn']
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
        vcf = "{yid}/%s" % config['callvnt2']['of35'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of35'],
    output:
        vcf = "{yid}/%s" % config['callvnt2']['of37'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of37'],
        stat = "{yid}/%s" % config['callvnt2']['of37'].replace(".vcf.gz", ".txt")
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        N = "%s.{yid}" % config['gatk_variant_filtration']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk_variant_filtration']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk_variant_filtration']['id']),
        j = lambda w: get_resource(w, config, 'gatk_variant_filtration'),
        mem = lambda w: get_resource(w, config, 'gatk_variant_filtration')['mem'] - 1,
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_variant_filtration')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        bcftools filter -i 'FILTER=="PASS"' {input.vcf} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        bcftools stats -s - {output.vcf} > {output.stat}
        """

rule gatk_snpeff:
    input:
        vcf = "{yid}/%s" % config['callvnt2']['of37'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of37'],
    output:
        vcf = "%s/{yid}/%s" % (config['oid'], config['callvnt2']['of38']),
        tbi = "%s/{yid}/%s.tbi" % (config['oid'], config['callvnt2']['of38']),
        stat = "%s/{yid}/%s" % (config['oid'], config['callvnt2']['of38'].replace(".vcf.gz", ".stat.txt"))
    params:
        genome = lambda w: config['y'][w.yid]['ref'],
        idx = lambda w: config['g'][config['y'][w.yid]['ref']]['snpeff']['xcfg'],
        stat0 = "{yid}/%s" % config['callvnt2']['of37'].replace(".vcf.gz", ".txt"),
        N = "%s.{yid}" % config['gatk_snpeff']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['gatk_snpeff']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['gatk_snpeff']['id']),
        j = lambda w: get_resource(w, config, 'gatk_snpeff'),
        mem = lambda w: get_resource(w, config, 'gatk_snpeff')['mem'] - 1,
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_snpeff')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.idx} {params.genome} -ud 0 \
                {input.vcf} | bgzip > {output.vcf}
        bcftools index -t {output.vcf}
        cp {params.stat0} {output.stat}
        """


