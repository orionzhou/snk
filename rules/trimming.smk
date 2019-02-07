rule fastp:
    input:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['trimming']['idir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['trimming']['idir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['trimming']['idir'],
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['fastp']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['fastp']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['fastp']['odir'],
        json = "{yid}/%s/{sid}.json" % config['fastp']['odir'],
        html = "{yid}/%s/{sid}.html" % config['fastp']['odir'],
    params:
        paired = lambda w: config['y'][w.yid]['t'][w.sid]['paired'],
        N = lambda w: "%s.%s.%s" % (w.yid, config['fastp']['id'], w.sid),
        e = lambda w: "%s/%s/%s/%s.e" % (w.yid, config['dirp'], config['fastp']['id'], w.sid),
        o = lambda w: "%s/%s/%s/%s.o" % (w.yid, config['dirp'], config['fastp']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fastp')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fastp')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fastp')['mem']
    threads: config['fastp']['ppn']
    run:
        if params.paired:
            shell("""
            fastp --thread {threads} \
            -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html}
            touch {output.r0}
            """)
        else:
            shell("""
            fastp --thread {threads} \
            -i {input.r0} -o {output.r0} \
            -j {output.json} -h {output.html}
            touch {output.r1} {output.r2}
            """)

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
        paired = lambda w: config['y'][w.yid]['t'][w.sid]['paired'],
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_pe'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"],
        N = lambda w: "%s.%s.%s" % (w.yid, config['trimmomatic']['id'], w.sid),
        e = lambda w: "%s/%s/%s/%s.e" % (w.yid, config['dirp'], config['trimmomatic']['id'], w.sid),
        o = lambda w: "%s/%s/%s/%s.o" % (w.yid, config['dirp'], config['trimmomatic']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['mem']
    threads: config['trimmomatic']['ppn']
    run:
        if params.paired:
            shell("""
            trimmomatic SE -threads {threads} \
            {input} {output} \
            {params.trimmer} >{log} 2>&1
            touch {output.r1} {output.r2} {output.r1u} {output.r2u}
            """)
        else:
            shell("""
            trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} {output.r1} {output.r1u} {output.r2} {output.r2u} \
            {params.trimmer} >{log} 2>&1
            touch {output.r0}
            """)

rule bbduk:
    input:
        "{yid}/%s/{sid}.fq.gz" % config['trimming']['idir']
    output:
        "{yid}/%s/{sid}.fq.gz" % config['bbduk']['odir'],
        "{yid}/%s/{sid}.se.log" % config['bbduk']['odir'],
    params:
        cmd = config['bbduk']['cmd'],
        paired = lambda w: config['y'][w.yid]['t'][w.sid]['paired'],
        extra = "ref=%s %s" %
            (','.join(config['bbduk']['refs']), config['bbduk']['extra']),
        N = lambda w: "%s.%s.%s" % (w.yid, config['bbduk']['id'], w.sid),
        e = lambda w: "%s/%s/%s/%s.e" % (w.yid, config['dirp'], config['bbduk']['id'], w.sid),
        o = lambda w: "%s/%s/%s/%s.o" % (w.yid, config['dirp'], config['bbduk']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'bbduk')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'bbduk')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'bbduk')['mem']
    threads: config['bbduk']['ppn']
    shell:
        "{params.cmd} in={input} out={output[0]} {params.extra} stats={output[1]}"

def trimming_inputs(w):
    yid, sid = w.yid, w.sid
    readtype = config['y'][yid]['readtype']
    assert readtype in config['valid']['readtype'], "invalid readtype: %s" % readtype
    idir = ''
    if readtype == '3rnaseq':
        idir = config['bbduk']['odir']
    elif readtype in ['illumina','solexa']:
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
    shell:
        "jsonutil.py fastp {input} > {output}"


