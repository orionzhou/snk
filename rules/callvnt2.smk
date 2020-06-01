include: 'gatk.smk'

#def genomics_dbimport_inputs(w):
def cv21_inputs(w):
    yid = w.yid
    vcfs, tbis = [], []
    for sid, sdic in config['y'][yid]['t'].items():
        oyid, gt = sdic['yid'], sdic['gt']
        fv = '%s/%s/%s.g.vcf.gz' % (config['callvnt2']['vdir'], oyid, gt)
        fx = '%s/%s/%s.g.vcf.gz.tbi' % (config['callvnt2']['vdir'], oyid, gt)
        vcfs.append(fv)
        tbis.append(fx)
    return {'vcfs':vcfs, 'tbis':tbis}

rule cv21_combine_gvcfs:
    input: unpack(cv21_inputs)
    output:
        vcf = "{yid}/%s/{rid}.g.vcf.gz" % config['callvnt2']['od05'],
        tbi = "{yid}/%s/{rid}.g.vcf.gz.tbi" % config['callvnt2']['od05'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win127'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['cv21_combine_gvcfs']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['cv21_combine_gvcfs']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['cv21_combine_gvcfs']['id']),
        j = lambda w: get_resource(w, config, 'cv21_combine_gvcfs'),
        mem = lambda w: get_resource(w, config, 'cv21_combine_gvcfs')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv21_combine_gvcfs')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" CombineGVCFs \
            {params.extra} -R {params.ref} -L {params.region} \
            {params.gvcfs} -O {output.vcf}
        """

rule cv22_genotype_gvcf:
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
        N = "{yid}.%s.{rid}" % config['cv22_genotype_gvcf']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['cv22_genotype_gvcf']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['cv22_genotype_gvcf']['id']),
        j = lambda w: get_resource(w, config, 'cv22_genotype_gvcf'),
        mem = lambda w: get_resource(w, config, 'cv22_genotype_gvcf')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv22_genotype_gvcf')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
#        -L {input.bed} --include-non-variant-sites \
#        -V gendb://{params.gendb} \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" GenotypeGVCFs \
            {params.extra} -R {params.ref} \
            --allow-old-rms-mapping-quality-annotation-data \
            -V {input.vcf} \
            -O {output.vcf}
        """

rule cv23_merge_vcfs:
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
        N = "{yid}.%s" % config['cv23_merge_vcfs']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv23_merge_vcfs']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv23_merge_vcfs']['id']),
        j = lambda w: get_resource(w, config, 'cv23_merge_vcfs'),
        mem = lambda w: get_resource(w, config, 'cv23_merge_vcfs')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv23_merge_vcfs')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
#        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
#            {params.extra} {params.input_str} -O {output.vcf}
        """
        bcftools concat {input.vcfs} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """


#___ variant recal
rule cv24a_pick_training:
    input:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt2']['od20'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt2']['od20'],
    output:
        snp = "{yid}/%s/{rid}.snp.vcf.gz" % config['callvnt2']['od31'],
        snp_tbi = "{yid}/%s/{rid}.snp.vcf.gz.tbi" % config['callvnt2']['od31'],
        idl = "{yid}/%s/{rid}.idl.vcf.gz" % config['callvnt2']['od31'],
        idl_tbi = "{yid}/%s/{rid}.idl.vcf.gz.tbi" % config['callvnt2']['od31']
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        minqual = 990,
        N = "{yid}.%s.{rid}" % config['cv24a_pick_training']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['cv24a_pick_training']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['cv24a_pick_training']['id']),
        j = lambda w: get_resource(w, config, 'cv24a_pick_training'),
        mem = lambda w: get_resource(w, config, 'cv24a_pick_training')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv24a_pick_training')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        bcftools view -G -i 'TYPE="snp" & QUAL>{params.minqual}' {input.vcf} -Oz -o {output.snp}
        bcftools view -G -i 'TYPE="indel" & QUAL>{params.minqual}' {input.vcf} -Oz -o {output.idl}
        bcftools index -t {output.snp}
        bcftools index -t {output.idl}
        """

rule cv24b_merge_vcfs:
    input:
        snps = lambda w: expand("%s/%s/{rid}.snp.vcf.gz" % (w.yid, config['callvnt2']['od31']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win127'].keys())),
        snp_tbis = lambda w: expand("%s/%s/{rid}.snp.vcf.gz.tbi" % (w.yid, config['callvnt2']['od31']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win127'].keys())),
        idls = lambda w: expand("%s/%s/{rid}.idl.vcf.gz" % (w.yid, config['callvnt2']['od31']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win127'].keys())),
        idl_tbis = lambda w: expand("%s/%s/{rid}.idl.vcf.gz.tbi" % (w.yid, config['callvnt2']['od31']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win127'].keys())),
    output:
        snp = "{yid}/%s" % config['callvnt2']['of31a'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31a'],
        idl = "{yid}/%s" % config['callvnt2']['of31b'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31b'],
    params:
        cmd = config['gatk']['cmd'],
        snps = lambda w, input: ["-I %s" % x for x in input.snps],
        idls = lambda w, input: ["-I %s" % x for x in input.idls],
        extra = gatk_extra(picard = True, jdk = False),
        N = "{yid}.%s" % config['cv24b_merge_vcfs']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv24b_merge_vcfs']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv24b_merge_vcfs']['id']),
        j = lambda w: get_resource(w, config, 'cv24b_merge_vcfs'),
        mem = lambda w: get_resource(w, config, 'cv24b_merge_vcfs')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv24b_merge_vcfs')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
#        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
#            {params.extra} {params.snps} -O {output.snp}
#        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
#            {params.extra} {params.idls} -O {output.idl}
        """
        bcftools concat {input.snps} -Oz -o {output.snp}
        bcftools index -t {output.snp}
        bcftools concat {input.idls} -Oz -o {output.idl}
        bcftools index -t {output.idl}
        bcftools stats {output.snp} > {output.snp}.txt
        bcftools stats {output.idl} > {output.idl}.txt
        """

rule cv25a_snp_recal:
    input:
        vcf = "{yid}/%s" % config['callvnt2']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of30'],
        snp = "{yid}/%s" % config['callvnt2']['of31a'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31a'],
    output:
        snp_recal = "{yid}/%s" % config['callvnt2']['of32a'],
        snp_tranch = "{yid}/%s" % config['callvnt2']['of32b'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "{yid}.%s" % config['cv25a_snp_recal']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv25a_snp_recal']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv25a_snp_recal']['id']),
        j = lambda w: get_resource(w, config, 'cv25a_snp_recal'),
        mem = lambda w: get_resource(w, config, 'cv25a_snp_recal')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv25a_snp_recal')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" VariantRecalibrator \
        {params.extra} -R {params.ref} \
        --resource:hmp3,known=false,training=true,truth=true,prior=10.0 {input.snp} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
        -mode SNP \
        -V {input.vcf} -O {output.snp_recal} --tranches-file {output.snp_tranch}
        """

rule cv25b_snp_vqsr:
    input:
        vcf = "{yid}/%s" % config['callvnt2']['of30'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of30'],
        snp_recal = "{yid}/%s" % config['callvnt2']['of32a'],
        snp_tranch = "{yid}/%s" % config['callvnt2']['of32b'],
    output:
        snp_vcf = "{yid}/%s" % config['callvnt2']['of33'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of33'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "{yid}.%s" % config['cv25b_snp_vqsr']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv25b_snp_vqsr']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv25b_snp_vqsr']['id']),
        j = lambda w: get_resource(w, config, 'cv25b_snp_vqsr'),
        mem = lambda w: get_resource(w, config, 'cv25b_snp_vqsr')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv25b_snp_vqsr')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyVQSR \
        {params.extra} -R {params.ref} \
        -V {input.vcf} -O {output.snp_vcf} \
        --recal-file {input.snp_recal} --tranches-file {input.snp_tranch} \
        --truth-sensitivity-filter-level 99.5 \
        -mode SNP \
        --create-output-variant-index true
        """

rule cv25c_idl_recal:
    input:
        snp_vcf = "{yid}/%s" % config['callvnt2']['of33'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of33'],
        idl = "{yid}/%s" % config['callvnt2']['of31b'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt2']['of31b'],
    output:
        idl_recal = "{yid}/%s" % config['callvnt2']['of34a'],
        idl_tranch = "{yid}/%s" % config['callvnt2']['of34b'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "{yid}.%s" % config['cv25c_idl_recal']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv25c_idl_recal']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv25c_idl_recal']['id']),
        j = lambda w: get_resource(w, config, 'cv25c_idl_recal'),
        mem = lambda w: get_resource(w, config, 'cv25c_idl_recal')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv25c_idl_recal')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" VariantRecalibrator \
        {params.extra} -R {params.ref} \
        --resource:hmp3,known=false,training=true,truth=true,prior=10.0 {input.idl} \
        -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
        -mode INDEL \
        -V {input.snp_vcf} -O {output.idl_recal} --tranches-file {output.idl_tranch}
        """

rule cv25d_idl_vqsr:
    input:
        snp_vcf = "{yid}/%s" % config['callvnt2']['of33'],
        snp_tbi = "{yid}/%s.tbi" % config['callvnt2']['of33'],
        idl_recal = "{yid}/%s" % config['callvnt2']['of34a'],
        idl_tranch = "{yid}/%s" % config['callvnt2']['of34b'],
    output:
        idl_vcf = "{yid}/%s" % config['callvnt2']['of35'],
        idl_tbi = "{yid}/%s.tbi" % config['callvnt2']['of35']
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(jdk = True),
        N = "{yid}.%s" % config['cv25d_idl_vqsr']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv25d_idl_vqsr']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv25d_idl_vqsr']['id']),
        j = lambda w: get_resource(w, config, 'cv25d_idl_vqsr'),
        mem = lambda w: get_resource(w, config, 'cv25d_idl_vqsr')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv25d_idl_vqsr')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyVQSR \
        {params.extra} -R {params.ref} \
        -V {input.snp_vcf} -O {output.idl_vcf} \
        --recal-file {input.idl_recal} --tranches-file {input.idl_tranch} \
        --truth-sensitivity-filter-level 99.0 \
        -mode INDEL \
        --create-output-variant-index true
        """

# final filter
rule cv27a_filter:
    input:
        vcf = "{yid}/%s" % config['callvnt2']['of35'],
        tbi = "{yid}/%s.tbi" % config['callvnt2']['of35'],
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt2']['od37'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt2']['od37'],
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['cv27a_filter']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['cv27a_filter']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['cv27a_filter']['id']),
        j = lambda w: get_resource(w, config, 'cv27a_filter'),
        mem = lambda w: get_resource(w, config, 'cv27a_filter')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv27a_filter')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
#bcftools filter -i 'FILTER=="PASS"' {input.vcf} -Oz -o {output.vcf}
#bcftools index -t {output.vcf}
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" SelectVariants \
            {params.extra} -R {params.ref} \
            -L {params.region} --exclude-filtered \
            -V {input.vcf} -O {output.vcf}
        """

rule cv28a_snpeff:
    input:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt2']['od37'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt2']['od37'],
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt2']['od38'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt2']['od38'],
    params:
        genome = lambda w: config['y'][w.yid]['ref'],
        idx = lambda w: config['g'][config['y'][w.yid]['ref']]['snpeff']['xcfg'],
        N = "{yid}.%s.{rid}" % config['cv28a_snpeff']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['cv28a_snpeff']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['cv28a_snpeff']['id']),
        j = lambda w: get_resource(w, config, 'cv28a_snpeff'),
        mem = lambda w: get_resource(w, config, 'cv28a_snpeff')['mem'] - 1,
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv28a_snpeff')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.idx} {params.genome} -ud 0 \
                {input.vcf} | bgzip > {output.vcf}
        bcftools index -t {output.vcf}
        """

rule cv28b_merge_vcfs:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.vcf.gz" % (w.yid, config['callvnt2']['od38']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
        tbis = lambda w: expand("%s/%s/{rid}.vcf.gz.tbi" % (w.yid, config['callvnt2']['od38']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
    output:
        vcf = "%s/{yid}.vcf.gz" % config['callvnt2']['odir'],
        tbi = "%s/{yid}.vcf.gz.tbi" % config['callvnt2']['odir'],
    params:
        cmd = config['gatk']['cmd'],
        genome = lambda w: config['y'][w.yid]['ref'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        vcfs = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True, jdk = False),
        N = "{yid}.%s" % config['cv28b_merge_vcfs']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv28b_merge_vcfs']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv28b_merge_vcfs']['id']),
        j = lambda w: get_resource(w, config, 'cv28b_merge_vcfs'),
        mem = lambda w: get_resource(w, config, 'cv28b_merge_vcfs')['mem'] - 1,
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv28b_merge_vcfs')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
#        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
#            {params.extra} {params.vcfs} -O {output.vcf}
        """
        bcftools concat {input.vcfs} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule cv29a_stat:
    input:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt2']['od37'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['callvnt2']['od37'],
    output:
        "{yid}/%s/{rid}.txt" % config['callvnt2']['od39'],
    params:
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        N = "{yid}.%s.{rid}" % config['cv29a_stat']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['cv29a_stat']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['cv29a_stat']['id']),
        j = lambda w: get_resource(w, config, 'cv29a_stat'),
        mem = lambda w: get_resource(w, config, 'cv29a_stat')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv29a_stat')['ppn']
    conda: "../envs/callvnt.yml"
    shell: "bcftools stats -s - {input.vcf} > {output}"

rule cv29b_merge_stats:
    input:
        lambda w: expand("%s/%s/{rid}.txt" % (w.yid, config['callvnt2']['od39']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
    output:
        "%s/{yid}/%s" % (config['oid'], config['callvnt2']['out_stat'])
    params:
        N = "{yid}.%s" % config['cv29b_merge_stats']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['cv29b_merge_stats']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['cv29b_merge_stats']['id']),
        j = lambda w: get_resource(w, config, 'cv29b_merge_stats'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv29b_merge_stats')['ppn']
    conda: "../envs/callvnt.yml"
    shell: 'merge.bcftools.stats.R {input} {output}'


