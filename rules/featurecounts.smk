def featurecounts_extra(wildcards):
    extra = config["featurecounts"]["extra"]
    extra += " --tmpDir %s" % config['tmpdir']
    extra += " -p"
    if config['stranded'] == 'reverse':
        extra += " -s 2"
    elif config['stranded'] == 'yes':
        extra += " -s 1"
    return extra

rule featurecounts:
    input:
        expand(["%s/{sid}.bam" % config['featurecounts']['idir']], sid = config['SampleID'])
    output:
        protected("%s/01.txt" % config['featurecounts']['odir']),
        protected("%s/01.txt.summary" % config['featurecounts']['odir'])
        
    log:
        "%s/featurecounts.log" % config['dirl']
    params:
        cmd = config['featurecounts']['cmd'],
        gtf = config['featurecounts']['gtf'],
        extra = featurecounts_extra,
    threads:
        config["featurecounts"]["threads"]
    shell:
        "{params.cmd} -T {threads} {params.extra} "
        "-a {params.gtf} -o {output[0]} {input} "
        ">{log} 2>&1"

