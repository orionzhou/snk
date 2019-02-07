def fq_compress_inputs(w):
    yid, sid = w.yid, w.sid
    inputs = dict()
    if config['y'][yid]['t'][sid]['paired']:
        inputs['r1'] = config['y'][yid]['t'][sid]['r1']
        inputs['r2'] = config['y'][yid]['t'][sid]['r2']
    else:
        inputs['r0'] = config['y'][yid]['t'][sid]['r0']
    return inputs

rule fq_compress:
    input: unpack(fq_compress_inputs)
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['fq_compress']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['fq_compress']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['fq_compress']['odir']
    params:
        paired = lambda w: config['y'][w.yid]['t'][w.sid]['paired'],
        N = lambda w: "%s.%s.%s" % (w.yid, config['fq_compress']['id'], w.sid),
        e = lambda w: "%s/%s/%s/%s.e" % (w.yid, config['dirp'], config['fq_compress']['id'], w.sid),
        o = lambda w: "%s/%s/%s/%s.o" % (w.yid, config['dirp'], config['fq_compress']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['mem']
    threads: config['fq_compress']['ppn']
    run:
        if params.paired:
            if input.r1.endswith(".gz"):
                shell("""
                ln -sf {input.r1} {output.r1}
                ln -sf {input.r2} {output.r2}
                touch {output.r0}
                """)
            else:
                shell("""
                pigz -p {threads} -c {input.r1} > {output.r1}
                pigz -p {threads} -c {input.r2} > {output.r2}
                touch {output.r0}
                """)
        else:
            if input.r0.endswith(".gz"):
                shell("""
                ln -sf {input.r0} {output.r0}
                touch {output.r1} {output.r2}
                """)
            else:
                shell("""
                pigz -p {threads} -c {input.r0} > {output.r0}
                touch {output.r1} {output.r2}
                """)

rule fq_deinterleave:
    input: lambda w: config['y'][w.yid]['t'][w.sid]['r0']
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['fq_deinterleave']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['fq_deinterleave']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['fq_deinterleave']['odir']
    params:
        extra = config["fq_deinterleave"]["extra"],
        N = lambda w: "{yid}.%s.{sid}" % config['fq_deinterleave']['id'],
        e = lambda w: "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['fq_deinterleave']['id']),
        o = lambda w: "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['fq_deinterleave']['id']),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fq_deinterleave')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fq_deinterleave')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fq_deinterleave')['mem']
    threads: config['fq_deinterleave']['ppn']
    shell:
        """
        zcat {input} | \
        deinterleave_fastq.sh {output.r1} {output.r2} {threads} compress
        touch {output.r0}
        """

rule fq_dump:
    input:
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['fq_dump']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['fq_dump']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['fq_dump']['odir']
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['fq_dump']['id'])
    params:
        paired = lambda w: config['y'][w.yid]['t'][w.sid]['paired'],
        odir = "{yid}/%s" % config['fq_dump']['odir'],
        o0 = "{yid}/%s/{sid}.fastq" % config['fq_dump']['odir'],
        o1 = "{yid}/%s/{sid}_1.fastq" % config['fq_dump']['odir'],
        o2 = "{yid}/%s/{sid}_2.fastq" % config['fq_dump']['odir'],
        tmp = config['tmpdir'],
        N = "{yid}.%s.{sid}" % config['fq_dump']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['fq_dump']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['fq_dump']['id']),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fq_dump')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fq_dump')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fq_dump')['mem'],
        load = lambda w, attempt:  get_resource(config, attempt, 'fq_dump')['load']
    threads: config['fq_dump']['ppn']
    run:
        if params.paired:
            shell("""
            fasterq-dump --split-files -e {threads} -m {params.mem} \
                    -O {params.odir} -t {params.tmp} {wildcards.sid} \
                    >{log} 2>&1

            pigz -p {threads} --fast -c {params.o1} >{output.r1}
            pigz -p {threads} --fast -c {params.o2} >{output.r2}
            rm {params.o1} {params.o2}
            touch {output.r0}
            """)
        else:
            shell("""
            fasterq-dump --split-files -e {threads} -m {params.mem} \
                    -O {params.odir} -t {params.tmp} {wildcards.sid} \
                    >{log} 2>&1

            pigz -p {threads} --fast -c {params.o0} >{output.r0}
            rm {params.o0}
            touch {output.r1} {output.r2}
            """)

def fastq_inputs(w):
    yid, sid = w.yid, w.sid
    source = config['y'][yid]['source']
    assert source in config['valid']['source'], "invalid source: %s" % source
    idir = ''
    if source == 'sra':
        idir = config['fq_dump']['odir']
    elif source == 'local':
        idir = config['fq_compress']['odir']
    elif source == 'local_interleaved':
        idir = config['fq_deinterleave']['odir']
    pre = op.abspath("%s/%s/%s" % (yid, idir, sid))
    pre = "%s/%s/%s" % (yid, idir, sid)
    inputs = {
        'r0': "%s.fq.gz" % pre,
        'r1': "%s_1.fq.gz" % pre,
        'r2': "%s_2.fq.gz" % pre,
    }
    return inputs

rule fastq:
    input: unpack(fastq_inputs)
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['fastq']['odir'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['fastq']['odir'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['fastq']['odir'],
    run:
        shell("ln -sf %s %s" % (op.abspath(input.r0), output.r0))
        shell("ln -sf %s %s" % (op.abspath(input.r1), output.r1))
        shell("ln -sf %s %s" % (op.abspath(input.r2), output.r2))


