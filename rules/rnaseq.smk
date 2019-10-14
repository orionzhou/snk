def featurecounts_extra(w):
    yid = w.yid
    extra = "--primary -Q 20 -t exon -g gene_id --byReadGroup"
    extra += " --tmpDir %s" % config['tmpdir']
    extra += " -p"
    if config['y'][yid]['stranded'] == 'reverse':
        extra += " -s 2"
    elif config['y'][yid]['stranded'] == 'yes':
        extra += " -s 1"
    return extra

rule featurecounts:
    input: "{yid}/%s/{sid}.bam" % config['rnaseq']['idir']
    output:
        "{yid}/%s/{sid}.txt" % config['rnaseq']['od31f'],
        "{yid}/%s/{sid}.txt.summary" % config['rnaseq']['od31f']
    params:
        gtf = lambda w: db_index(w, 'gtf'),
        extra = featurecounts_extra,
        N = "{yid}.%s.{sid}" % config['featurecounts']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['featurecounts']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['featurecounts']['id']),
        j = lambda w: get_resource(w, config, 'featurecounts'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'featurecounts')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        featureCounts -T {threads} {params.extra} \
            -a {params.gtf} -o {output[0]} {input}
        """

rule merge_featurecounts:
    input:
        lambda w: expand("%s/%s/{sid}.txt" % (w.yid, config['rnaseq']['od31f']), sid = config['y'][w.yid]['SampleID'])
    output: protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_fcnt']))
    params:
        N = "{yid}.%s" % (config['merge_featurecounts']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['merge_featurecounts']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['merge_featurecounts']['id']),
        j = lambda w: get_resource(w, config, 'merge_featurecounts'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'merge_featurecounts')['ppn']
    conda: "../envs/r.yml"
    shell: "merge.stats.R --opt featurecounts -o {output} {input}"

def mmquant_strand(w):
    yid = w.yid
    if config['y'][yid]['stranded'] == 'reverse':
        return "RF"
    elif config['y'][yid]['stranded'] == 'yes':
        return "FR"
    else:
        return "U"

rule mmquant:
    input: "{yid}/%s/{sid}.bam" % config['rnaseq']['idir']
    output:
        "{yid}/%s/{sid}.tsv" % config['rnaseq']['od31m'],
        "{yid}/%s/{sid}.txt" % config['rnaseq']['od31m']
    params:
        gtf = lambda w: db_index(w, 'gtf'),
#        gtf = '/home/springer/zhoux379/projects/genome/data/B73xMo17/15_intervals/02.chrom.gtf',
        strand = mmquant_strand,
        N = "{yid}.%s.{sid}" % config['mmquant']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['mmquant']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['mmquant']['id']),
        j = lambda w: get_resource(w, config, 'mmquant'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'mmquant')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        mmquant -a {params.gtf} -r {input} \
            -o {output[0]} -O {output[1]} \
            -l 0.3 -t {threads} -s {params.strand}
        """

rule merge_mmquant:
    input:
        lambda w: expand("%s/%s/{sid}.tsv" % (w.yid, config['rnaseq']['od31m']), sid = config['y'][w.yid]['SampleID'])
    output: protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_mmq']))
    params:
        N = "{yid}.%s" % (config['merge_mmquant']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['merge_mmquant']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['merge_mmquant']['id']),
        j = lambda w: get_resource(w, config, 'merge_mmquant'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'merge_mmquant')['ppn']
    conda: "../envs/r.yml"
    shell: "merge.stats.R --opt mmquant -o {output} {input}"

rule rc2cpm:
    input:
        sam = lambda w: config['y'][w.yid]['samplelist'],
        exp = "%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_fcnt']),
        cfg = lambda w: config['g'][config['y'][w.yid]['ref']]["rds"]['xout']
    output: protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_cpm']))
    params:
        N = "{yid}.%s" % (config['rc2cpm']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['rc2cpm']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['rc2cpm']['id']),
        j = lambda w: get_resource(w, config, 'rc2cpm'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'rc2cpm')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        rc2cpm.R {input.sam} {input.exp} {output} \
            --opt featurecounts --yid {wildcards.yid} --config {input.cfg}
        """

rule rc2cpm2:
    input:
        sam = "%s/11_qc/{yid}/%s" % (config['dirh'], config['rnaseq']['meta']),
        exp = "%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_fcnt']),
        cfg = lambda w: config['g'][config['y'][w.yid]['ref']]["rds"]['xout']
    output: protected("%s/11_qc/{yid}/%s" % (config['dirh'], config['rnaseq']['out_cpm']))
    params:
        N = "{yid}.%s" % (config['rc2cpm']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['rc2cpm']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['rc2cpm']['id']),
        j = lambda w: get_resource(w, config, 'rc2cpm'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'rc2cpm')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        rc2cpm.R {input.sam} {input.exp} {output} \
            --opt featurecounts --yid {wildcards.yid} --config {input.cfg}
        """


