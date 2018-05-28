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
        expand(["22.bam/{sid}_{gt}.bam"],
            zip, sid = config['t']['sid'], gt = config['t']['genotype'])
    output:
        "24.featurecounts/01.txt"
    log:
        "logs/featurecounts.log"
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
        "{input} ")
        #"> {log}")

