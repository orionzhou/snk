def featurecounts_extra(wildcards):
    extra = config["featurecounts"]["extra"]
    extra += " --tmpDir %s" % config['tmpdir']
    if config['paired']:
        extra += " -p"
    if config['stranded'] == 'reverse':
        extra += " -s 2"
    elif config['stranded'] == 'yes':
        extra += " -s 1"
    return extra

rule featurecounts:
    input:
        expand(["%s/{sid}_{gt}.bam" % config['featurecounts']['idir']],
            zip, sid = config['t']['SampleID'], gt = config['t']['Genotype'])
    output:
        config['featurecounts']['outfile']
    log:
        "%s/featurecounts.log" % config['dirl']
    params:
        cmd = config['featurecounts']['cmd'],
        gtf = config['featurecounts']['gtf'],
        extra = featurecounts_extra,
    threads:
        config["featurecounts"]["threads"]
    run:
        shell("{params.cmd} "
        "-T {threads} "
        "{params.extra} "
        "-a {params.gtf} "
        "-o {output} "
        "{input} "
        "> {log}")

