def featurecounts_extra(w):
    yid = w.yid
    extra = config["featurecounts"]["extra"]
    extra += " --tmpDir %s" % config['tmpdir']
    extra += " -p"
    if config['y'][yid]['stranded'] == 'reverse':
        extra += " -s 2"
    elif config['y'][yid]['stranded'] == 'yes':
        extra += " -s 1"
    return extra

rule featurecounts:
    input:
        "{yid}/%s/{sid}.bam" % config['featurecounts']['idir']
    output:
        "{yid}/%s/{sid}.txt" % config['featurecounts']['odir'],
        "{yid}/%s/{sid}.txt.summary" % config['featurecounts']['odir']
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['featurecounts']['id'])
    params:
        gtf = lambda w: config[config['y'][w.yid]['reference']]['gtf'],
        extra = featurecounts_extra,
        N = "{yid}.%s.{sid}" % config['featurecounts']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['featurecounts']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['featurecounts']['id'])
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'featurecounts')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'featurecounts')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'featurecounts')['mem']
    threads: config['featurecounts']['ppn']
    conda: "../envs/job.yml"
    shell:
        "featureCounts -T {threads} {params.extra} "
        "-a {params.gtf} -o {output[0]} {input} "
        ">{log} 2>&1"

def merge_featurecounts_inputs(w):
    yid = w.yid
    inputs = []
    for sid in config['y'][yid]['SampleID']:
        inputs.append("%s/%s/%s.txt" % (yid, config['featurecounts']['odir'], sid))
    return inputs

rule merge_featurecounts:
    input: merge_featurecounts_inputs
    output:
        protected("{yid}/%s/%s" % (config['dird'], config['merge_featurecounts']['out']))
    params:
        N = "{yid}.%s" % (config['merge_featurecounts']['id']),
        e = "{yid}/%s/%s.e" % (config['dirp'], config['merge_featurecounts']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['merge_featurecounts']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'merge_featurecounts')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'merge_featurecounts')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'merge_featurecounts')['mem']
    threads: config['merge_featurecounts']['ppn']
    conda: "../envs/r.yml"
    shell:
        "merge.featurecounts.R -o {output} {input}"

rule rc2cpm_raw:
    input:
        sam = lambda w: config['y'][w.yid]['samplelist'],
        exp = "{yid}/%s/%s" % (config['dird'], config['merge_featurecounts']['out']),
        cfg = lambda w: config[config['y'][w.yid]['reference']]["rds"],
    output:
        protected("{yid}/%s/%s" % (config['dird'], config['rc2cpm']['out_raw']))
    params:
        N = "{yid}.%s" % (config['rc2cpm']['id']),
        e = "{yid}/%s/%s.e" % (config['dirp'], config['rc2cpm']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['rc2cpm']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'rc2cpm')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'rc2cpm')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'rc2cpm')['mem']
    threads: config['rc2cpm']['ppn']
    conda: "../envs/r.yml"
    shell:
        "rc2cpm.R {input.sam} {input.exp} {output} --yid {wildcards.yid} --config {input.cfg}"

rule rc2cpm:
    input:
        sam = lambda w: config['y'][w.yid]['samplelistc'],
        exp = "{yid}/%s/%s" % (config['dird'], config['merge_featurecounts']['out']),
        cfg = lambda w: config[config['y'][w.yid]['reference']]["rds"],
    output:
        protected("{yid}/%s/%s" % (config['dird'], config['rc2cpm']['out']))
    params:
        N = "{yid}.%s" % (config['rc2cpm']['id']),
        e = "{yid}/%s/%s.e" % (config['dirp'], config['rc2cpm']['id']),
        o = "{yid}/%s/%s.o" % (config['dirp'], config['rc2cpm']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'rc2cpm')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'rc2cpm')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'rc2cpm')['mem']
    threads: config['rc2cpm']['ppn']
    conda: "../envs/r.yml"
    shell:
        "rc2cpm.R {input.sam} {input.exp} {output} --yid {wildcards.yid} --config {input.cfg}"


