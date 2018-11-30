def gatk_extra(picard = False, jdk = False, hc = False):
    extra = ''
    if picard:
        extra += ' --TMP_DIR %s' % config['tmpdir']
    else:
        extra += ' --tmp-dir %s' % config['tmpdir']
    if jdk:
        if picard:
            extra += ' --USE_JDK_DEFLATER --USE_JDK_INFLATER'
        else:
            extra += ' --use-jdk-deflater --use-jdk-inflater'
    if hc:
        extra += ' -pairHMM LOGLESS_CACHING'
    return extra

rule gatk_mark_duplicates:
    input:
        "%s/{sid}.bam" % config['cleanbam']['idir']
    output:
        "%s/{sid}.bam" % config['cleanbam']['odir1'],
        "%s/{sid}.dedup.txt" % config['cleanbam']['odir1'],
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['gatk']['mark_duplicates']['id'])
    params:
        cmd = config['gatk']['cmd'],
        extra = gatk_extra(picard = True, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['mark_duplicates']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['mark_duplicates']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['mark_duplicates']['id'], w.sid),
        q = lambda w, resources: resources.q,
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['mem']
    threads: config['gatk']['mark_duplicates']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MarkDuplicates \
        {params.extra} \
        -I {input} -O {output[0]} \
        -M {output[1]} \
        >>{log} 2>&1
        """

rule gatk_base_recalibrator:
    input:
        "%s/{sid}.bam" % config['cleanbam']['odir1']
    output:
        "%s/{sid}.table" % config['cleanbam']['odir2'],
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['gatk']['base_recalibrator']['id'])
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        vcf = config[config['reference']]['gatk']['known_sites'],
        extra = gatk_extra(picard = False, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['base_recalibrator']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['base_recalibrator']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['base_recalibrator']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['mem']
    threads: config['gatk']['base_recalibrator']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" BaseRecalibrator \
        {params.extra} \
        -R {params.ref} \
        -I {input} \
        --known-sites {params.vcf} \
        -O {output[0]} \
        >>{log} 2>&1
        """

rule gatk_apply_bqsr:
    input:
        "%s/{sid}.bam" % config['cleanbam']['odir1'],
        "%s/{sid}.table" % config['cleanbam']['odir2'],
    output:
        "%s/{sid}.bam" % config['cleanbam']['odir2'],
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['gatk']['apply_bqsr']['id'])
    params:
        cmd = config['gatk']['cmd'],
        ref = config[config['reference']]['gatk']['xref'],
        extra = gatk_extra(picard = False, jdk = True),
        N = lambda w: "%s.%s" % (config['gatk']['apply_bqsr']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk']['apply_bqsr']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk']['apply_bqsr']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['mem']
    threads: config['gatk']['apply_bqsr']['ppn']
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyBQSR \
        {params.extra} \
        -R {params.ref} \
        -I {input[0]} \
        --bqsr-recal-file {input[1]} \
        -O {output[0]} \
        >>{log} 2>&1
        """


