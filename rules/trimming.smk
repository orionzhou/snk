rule fastp:
    input:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['trimming']['idir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['trimming']['idir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['trimming']['idir']
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['fastp']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['fastp']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['fastp']['odir'],
        json = "{yid}/%s/{sid}.json" % config['fastp']['odir'],
        html = "{yid}/%s/{sid}.html" % config['fastp']['odir']
    params:
        paired = lambda w: int(config['y'][w.yid]['t'][w.sid]['paired']),
        N = "{yid}.%s.{sid}" % config['fastp']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['fastp']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['fastp']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fastp')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fastp')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fastp')['mem']
    threads: config['fastp']['ppn']
    conda: "../envs/job.yml"
    script: "../scripts/fastp.py"

rule trimmomatic:
    input:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['trimming']['idir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['trimming']['idir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['trimming']['idir'],
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['trimmomatic']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['trimmomatic']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['trimmomatic']['odir'],
        r1u = "{yid}/%s/{sid}_1.unpaired.fq.gz" % config['trimmomatic']['odir'],
        r2u = "{yid}/%s/{sid}_2.unpaired.fq.gz" % config['trimmomatic']['odir']
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['trimmomatic']['id'])
    params:
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_pe'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"],
        N = "{yid}.%s.{sid}" % config['trimmomatic']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['trimmomatic']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['trimmomatic']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['mem']
    threads: config['trimmomatic']['ppn']
    conda: "../envs/job.yml"
    script: "../scripts/trimmomatic.py"

rule bbduk:
    input:
        "{yid}/%s/{sid}.fq.gz" % config['trimming']['idir']
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['bbduk']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['bbduk']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['bbduk']['odir'],
        json = "{yid}/%s/{sid}.json" % config['bbduk']['odir'],
    params:
        cmd = config['bbduk']['cmd'],
        extra = "ref=%s %s" %
            (','.join(config['bbduk']['refs']), config['bbduk']['extra']),
        N = "{yid}.%s.{sid}" % config['bbduk']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['bbduk']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['bbduk']['id'])
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'bbduk')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'bbduk')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'bbduk')['mem']
    threads: config['bbduk']['ppn']
    conda: "../envs/job.yml"
    shell:
        "{params.cmd} in={input} out={output.r0} {params.extra} stats={output.json}"

def trimming_inputs(w):
    yid, sid = w.yid, w.sid
    readtype = config['y'][yid]['readtype']
    assert readtype in config['valid']['readtype'], "invalid readtype: %s" % readtype
    idir = ''
    if readtype == '3rnaseq':
        idir = config['bbduk']['odir']
        idir = config['fastp']['odir']
    elif readtype in ['illumina','solid']:
        idir = config['fastp']['odir']
    pre = op.abspath("%s/%s/%s" % (yid, idir, sid))
    pre = "%s/%s/%s" % (yid, idir, sid)
    inputs = {
        'r0': "%s.fq.gz" % pre,
        'r1': "%s_1.fq.gz" % pre,
        'r2': "%s_2.fq.gz" % pre,
        'json': "%s.json" % pre
    }
    return inputs

rule trimming:
    input: unpack(trimming_inputs)
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['trimming']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['trimming']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['trimming']['odir'],
        json = "{yid}/%s/{sid}.json" % config['trimming']['odir'],
    run:
        shell("ln -sf %s %s" % (op.abspath(input.r0), output.r0))
        shell("ln -sf %s %s" % (op.abspath(input.r1), output.r1))
        shell("ln -sf %s %s" % (op.abspath(input.r2), output.r2))
        shell("ln -sf %s %s" % (op.abspath(input.json), output.json))

def merge_trimstats_inputs(w):
    yid = w.yid
    inputs = []
    for sid in config['y'][yid]['SampleID']:
        pair_suf = 'pe' if config['y'][yid]['t'][sid]['paired'] else 'se'
        inputs.append("%s/%s/%s.json" % (yid, config['trimming']['odir'], sid))
    return inputs

rule merge_trimstats:
    input: merge_trimstats_inputs
    output:
        protected("{yid}/%s/%s" % (config['dird'], config['merge_trimstats']['out']))
    conda: "../envs/python.yml"
    shell:
        "jsonutil.py fastp {input} > {output}"


