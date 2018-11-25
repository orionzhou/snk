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
        "%s/%s/{sid}.log" % (config['dirl'], config['featurecounts']['id'])
    params:
        gtf = config[config['reference']]['gtf'],
        extra = featurecounts_extra,
        N = lambda w: "%s.%s" % (config['featurecounts']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['featurecounts']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['featurecounts']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'featurecounts')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'featurecounts')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'featurecounts')['mem']
    threads: config['featurecounts']['ppn']
    shell:
        "featureCounts -T {threads} {params.extra} "
        "-a {params.gtf} -o {output[0]} {input} "
        ">{log} 2>&1"

