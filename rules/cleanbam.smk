include: 'gatk.smk'

rule bam6_markdup:
    input: ancient("{yid}/%s/{sid}.bam" % config['cleanbam']['idir'])
    output:
        protected("{yid}/%s/{sid}.bam" % config['cleanbam']['od23']),
        protected("{yid}/%s/{sid}.dedup.txt" % config['cleanbam']['od23']),
    params:
        cmd = config['gatk']['cmd'],
        tmp = config['tmpdir'],
        extra = gatk_extra(picard = True, jdk = True),
        N = "{yid}.%s.{sid}" % config['bam6_markdup']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bam6_markdup']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bam6_markdup']['id']),
        j = lambda w: get_resource(w, config, 'bam6_markdup'),
        mem = lambda w: get_resource(w, config, 'bam6_markdup')['mem'],
    resources:
        attempt = lambda w, attempt: attempt,
        load = lambda w: get_resource(w, config, 'bam6_markdup')['load'],
    threads: lambda w: get_resource(w, config, 'bam6_markdup')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        picard -Xmx{params.mem}G MarkDuplicates \
        TMP_DIR={params.tmp} USE_JDK_DEFLATER=true USE_JDK_INFLATER=true \
        I={input} O={output[0]} M={output[1]}
        samtools index {output[0]}
        """
#        {params.cmd} --java-options "-Xmx{params.mem}G" MarkDuplicates \
#        {params.extra} \
#        -I {input} -O {output[0]} -M {output[1]} \
#        >>{log} 2>&1

def gatk_known_sites(w):
    known_sites = "%s/%s" % (w.yid, config['cleanbam']['of24c'])
    if 'known_sites' in config['g'][config['y'][w.yid]['ref']]['gatk']:
        known_sites = config['g'][config['y'][w.yid]['ref']]['gatk']['known_sites']
    return known_sites

rule bam7a_recal:
    input:
        bam = "{yid}/%s/{sid}.bam" % config['cleanbam']['od23'],
        known_sites = gatk_known_sites,
    output:
        protected("{yid}/%s/{sid}.table" % config['cleanbam']['od24']),
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{sid}" % config['bam7a_recal']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bam7a_recal']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bam7a_recal']['id']),
        j = lambda w: get_resource(w, config, 'bam7a_recal'),
        mem = lambda w: get_resource(w, config, 'bam7a_recal')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bam7a_recal')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" BaseRecalibrator \
        {params.extra} --known-sites {input.known_sites} \
        -R {params.ref} -I {input.bam} -O {output[0]}
        """

rule bam7b_bqsr:
    input:
        "{yid}/%s/{sid}.bam" % config['cleanbam']['od23'],
        "{yid}/%s/{sid}.table" % config['cleanbam']['od24'],
    output:
        protected("{yid}/%s/{sid}.bam" % config['cleanbam']['od24']),
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{sid}" % config['bam7b_bqsr']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bam7b_bqsr']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bam7b_bqsr']['id']),
        j = lambda w: get_resource(w, config, 'bam7b_bqsr'),
        mem = lambda w: get_resource(w, config, 'bam7b_bqsr')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bam7b_bqsr')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" ApplyBQSR \
        {params.extra} --bqsr-recal-file {input[1]} \
        -R {params.ref} -I {input[0]} -O {output[0]}
        """

rule bam8_stat:
    input: "{yid}/%s/{sid}.bam" % config['cleanbam']['od24']
    output: "{yid}/%s/{sid}.tsv" % config['cleanbam']['od24']
    params:
        N = "{yid}.%s.{sid}" % config['bam8_stat']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bam8_stat']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bam8_stat']['id']),
        j = lambda w: get_resource(w, config, 'bam8_stat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bam8_stat')['ppn']
    conda: "../envs/work.yml"
    shell: "bam.py stat {input} > {output}"

rule bam9_merge_stats:
    input:
        lambda w: expand(ancient("%s/%s/{sid}.tsv" % (w.yid, config['cleanbam']['od24'])), sid = config['y'][w.yid]['t'].keys())
    output:
        protected("%s/{yid}/%s" % (config['oid'], config['merge_bamstats']['outv']))
    params:
        N = "{yid}.%s" % (config['bam9_merge_stats']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['bam9_merge_stats']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['bam9_merge_stats']['id']),
        j = lambda w: get_resource(w, config, 'bam9_merge_stats'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bam9_merge_stats')['ppn']
    conda: "../envs/work.yml"
    shell: "merge.stats.R --opt bam_stat -o {output} {input}"



