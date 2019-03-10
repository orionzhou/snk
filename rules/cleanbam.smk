def gatk_extra(picard = False, jdk = False, hc = False):
    extra = ''
    if picard:
        extra += ' --TMP_DIR %s' % config['tmpdir'] #' TMP_DIR=%s' % config['tmpdir']
    else:
        extra += ' --tmp-dir %s' % config['tmpdir']
    if jdk:
        if picard:
            extra += ' --USE_JDK_DEFLATER --USE_JDK_INFLATER' #' USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
        else:
            extra += ' --use-jdk-deflater --use-jdk-inflater'
    if hc:
        extra += ' -pairHMM LOGLESS_CACHING'
    return extra

rule gatk_mark_duplicates:
    input:
        "{yid}/%s/{sid}.bam" % config['mapping']['odir']
    output:
        protected("{yid}/%s/{sid}.bam" % config['cleanbam']['odir1']),
        protected("{yid}/%s/{sid}.dedup.txt" % config['cleanbam']['odir1']),
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['gatk']['mark_duplicates']['id'])
    params:
        cmd = config['gatk']['cmd'],
        tmp = config['tmpdir'],
        extra = gatk_extra(picard = True, jdk = True),
        N = "{yid}.%s.{sid}" % config['gatk']['mark_duplicates']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['gatk']['mark_duplicates']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['gatk']['mark_duplicates']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['mem'] - 2,
        load = lambda w, attempt:  get_resource(config, attempt, 'gatk','mark_duplicates')['load']
    threads: config['gatk']['mark_duplicates']['ppn']
    conda: "../envs/gatk.yml"
    shell:
        """
        picard -Xmx{params.mem}G MarkDuplicates \
        TMP_DIR={params.tmp} USE_JDK_DEFLATER=true USE_JDK_INFLATER=true \
        I={input} O={output[0]} M={output[1]}
        """
#        {params.cmd} --java-options "-Xmx{params.mem}G" MarkDuplicates \
#        {params.extra} \
#        -I {input} -O {output[0]} -M {output[1]} \
#        >>{log} 2>&1

rule gatk_base_recalibrator:
    input:
        "{yid}/%s/{sid}.bam" % config['cleanbam']['odir1']
    output:
        protected("{yid}/%s/{sid}.table" % config['cleanbam']['odir2']),
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['gatk']['base_recalibrator']['id'])
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        vcf = lambda w: config[config['y'][w.yid]['reference']]['gatk']['known_sites'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{sid}" % config['gatk']['base_recalibrator']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['gatk']['base_recalibrator']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['gatk']['base_recalibrator']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['mem']
    threads: config['gatk']['base_recalibrator']['ppn']
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" BaseRecalibrator \
        {params.extra} --known-sites {params.vcf} \
        -R {params.ref} -I {input} -O {output[0]} \
        >>{log} 2>&1
        """

rule gatk_apply_bqsr:
    input:
        "{yid}/%s/{sid}.bam" % config['cleanbam']['odir1'],
        "{yid}/%s/{sid}.table" % config['cleanbam']['odir2'],
    output:
        protected("{yid}/%s/{sid}.bam" % config['cleanbam']['odir2']),
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['gatk']['apply_bqsr']['id'])
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config[config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{sid}" % config['gatk']['apply_bqsr']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['gatk']['apply_bqsr']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['gatk']['apply_bqsr']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['mem']
    threads: config['gatk']['apply_bqsr']['ppn']
    conda: "../envs/gatk.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyBQSR \
        {params.extra} --bqsr-recal-file {input[1]} \
        -R {params.ref} -I {input[0]} -O {output[0]} \
        >>{log} 2>&1
        """

rule bam_stat2:
    input:
        "{yid}/%s/{sid}.bam" % config['cleanbam']['odir2']
    output:
        "{yid}/%s/{sid}.tsv" % config['cleanbam']['odir2']
    params:
        N = "{yid}.%s.{sid}" % config['bam_stat']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['bam_stat']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['bam_stat']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['mem']
    threads: config['bam_stat']['ppn']
    conda: "../envs/python.yml"
    shell:
        "bam.py stat {input} > {output}"

rule merge_bamstats2:
    input:
        lambda w: expand("%s/%s/{sid}.tsv" % (w.yid, config['cleanbam']['odir2']), sid = config['y'][w.yid]['t'].keys())
    output:
        protected("{yid}/%s/%s" % (config['dird'], config['merge_bamstats']['outv']))
    conda: "../envs/r.yml"
    shell:
        "merge.bamstats.R -o {output} {input}"



