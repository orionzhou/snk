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
        "%s/{sid}.bam" % config['featurecounts']['idir']
    output:
        protected("%s/{sid}.txt" % config['featurecounts']['odir']),
        protected("%s/{sid}.txt.summary" % config['featurecounts']['odir'])
    log:
        "%s/featurecounts/{sid}.log" % config['dirl']
    params:
        cmd = config['featurecounts']['cmd'],
        gtf = config['genomes'][config['reference']]['gtf'],
        extra = featurecounts_extra,
        ppn = config['featurecounts']['ppn'],
        walltime = config['featurecounts']['walltime'],
        mem = config['featurecounts']['mem']
    threads: config['featurecounts']['ppn']
    shell:
        "{params.cmd} -T {threads} {params.extra} "
        "-a {params.gtf} -o {output[0]} {input} "
        ">{log} 2>&1"

rule merge_featurecounts:
    input:
        expand(["%s/{sid}.txt" % config['featurecounts']['odir']], sid = config['SampleID'])
    output:
        protected("%s/%s" % (config['dird'], config['featurecounts']['out']))
    shell:
        "merge.featurecounts.R -o {output} {input}"
