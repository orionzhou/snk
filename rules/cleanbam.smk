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
    input: "{yid}/%s/{sid}.bam" % config['mapping']['odir']
    output:
        protected("{yid}/%s/{sid}.bam" % config['cleanbam']['od23']),
        protected("{yid}/%s/{sid}.dedup.txt" % config['cleanbam']['od23']),
    params:
        cmd = config['gatk']['cmd'],
        tmp = config['tmpdir'],
        extra = gatk_extra(picard = True, jdk = True),
        N = "{yid}.%s.{sid}" % config['gatk']['mark_duplicates']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['gatk']['mark_duplicates']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['gatk']['mark_duplicates']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'mark_duplicates')['mem'] - 2,
        load = lambda w, attempt:  get_resource(config, attempt, 'gatk','mark_duplicates')['load']
    threads: config['gatk']['mark_duplicates']['ppn']
    conda: "../envs/work.yml"
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

def gatk_known_sites(w):
    known_sites = "%s/%s" % (w.yid, config['cleanbam']['of24c'])
    if 'known_sites' in config['g'][config['y'][w.yid]['reference']]['gatk']:
        known_sites = config['g'][config['y'][w.yid]['reference']]['gatk']['known_sites']
    return known_sites

rule gatk_base_recalibrator:
    input:
        bam = "{yid}/%s/{sid}.bam" % config['cleanbam']['od23'],
        known_sites = gatk_known_sites,
    output:
        protected("{yid}/%s/{sid}.table" % config['cleanbam']['od24']),
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{sid}" % config['gatk']['base_recalibrator']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['gatk']['base_recalibrator']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['gatk']['base_recalibrator']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'base_recalibrator')['mem']
    threads: config['gatk']['base_recalibrator']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" BaseRecalibrator \
        {params.extra} --known-sites {input.known_sites} \
        -R {params.ref} -I {input.bam} -O {output[0]}
        """

rule gatk_apply_bqsr:
    input:
        "{yid}/%s/{sid}.bam" % config['cleanbam']['od23'],
        "{yid}/%s/{sid}.table" % config['cleanbam']['od24'],
    output:
        protected("{yid}/%s/{sid}.bam" % config['cleanbam']['od24']),
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['gatk']['xref'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{sid}" % config['gatk']['apply_bqsr']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['gatk']['apply_bqsr']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['gatk']['apply_bqsr']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk', 'apply_bqsr')['mem']
    threads: config['gatk']['apply_bqsr']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyBQSR \
        {params.extra} --bqsr-recal-file {input[1]} \
        -R {params.ref} -I {input[0]} -O {output[0]}
        """

rule bam_stat2:
    input:
        "{yid}/%s/{sid}.bam" % config['cleanbam']['od24']
    output:
        "{yid}/%s/{sid}.tsv" % config['cleanbam']['od24']
    params:
        N = "{yid}.%s.{sid}" % config['bam_stat']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bam_stat']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bam_stat']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['mem']
    threads: config['bam_stat']['ppn']
    conda: "../envs/work.yml"
    shell:
        "bam.py stat {input} > {output}"

rule merge_bamstats2:
    input:
        lambda w: expand("%s/%s/{sid}.tsv" % (w.yid, config['cleanbam']['od24']), sid = config['y'][w.yid]['t'].keys())
    output:
        protected("{yid}/data/%s" % config['merge_bamstats']['outv'])
    conda: "../envs/work.yml"
    shell:
        "merge.bamstats.R -o {output} {input}"



