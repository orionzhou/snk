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
        r0 = "{yid}/%s/{sid}.fq.gz" % config['barn']['od09z'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['barn']['od09z'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['barn']['od09z']
    params:
        N = "{yid}.%s.{sid}" % config['fq_compress']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['fq_compress']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['fq_compress']['id']),
        j = lambda w: get_resource(w, config, 'fq_compress'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fq_compress')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/fq_compress.py"

rule fq_deinterleave:
    input: lambda w: config['y'][w.yid]['t'][w.sid]['r0']
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['barn']['od09v'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['barn']['od09v'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['barn']['od09v']
    params:
        N = "{yid}.%s.{sid}" % config['fq_deinterleave']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['fq_deinterleave']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['fq_deinterleave']['id']),
        j = lambda w: get_resource(w, config, 'fq_deinterleave'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fq_deinterleave')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        zcat {input} | \
            deinterleave_fastq.sh {output.r1} {output.r2} {threads} compress
        touch {output.r0}
        """

rule fq_dump:
    input:
    output:
        r0 = "{yid}/%s/{sid}.fq.gz" % config['barn']['od09d'],
        r1 = "{yid}/%s/{sid}_1.fq.gz" % config['barn']['od09d'],
        r2 = "{yid}/%s/{sid}_2.fq.gz" % config['barn']['od09d']
    params:
        odir = "{yid}/%s" % config['barn']['od09d'],
        o0 = "{yid}/%s/{sid}.fastq" % config['barn']['od09d'],
        o1 = "{yid}/%s/{sid}_1.fastq" % config['barn']['od09d'],
        o2 = "{yid}/%s/{sid}_2.fastq" % config['barn']['od09d'],
        tmp = config['tmpdir'],
        N = "{yid}.%s.{sid}" % config['fq_dump']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['fq_dump']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['fq_dump']['id']),
        j = lambda w: get_resource(w, config, 'fq_dump'),
        mem = lambda w: get_resource(w, config, 'fq_dump')['mem'],
    resources:
        attempt = lambda w, attempt: attempt,
        load = lambda w: get_resource(w, config, 'fq_dump')['load'],
    threads: lambda w: get_resource(w, config, 'fq_dump')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/fq_dump.py"

def fq_copy_inputs(w):
    yid, sid = w.yid, w.sid
    source = config['y'][yid]['source']
    assert source in config['valid']['source'], "invalid source: %s" % source
    idir = ''
    if source == 'sra':
        idir = config['barn']['od09d']
    elif source == 'local':
        fmt = config['y'][yid]['format']
        if config['y'][yid]['interleaved']:
            idir = config['barn']['od09v']
        else:
            idir = config['barn']['od09z']
    else:
        print("unsupported source: %s for %s" % (source, sid))
        sys.exit(1)
    pre = op.abspath("%s/%s/%s" % (yid, idir, sid))
    pre = "%s/%s/%s" % (yid, idir, sid)
    return dict(
        r0 = ancient("%s.fq.gz" % pre),
        r1 = ancient("%s_1.fq.gz" % pre),
        r2 = ancient("%s_2.fq.gz" % pre)
    )

rule fq_copy:
    input: unpack(fq_copy_inputs)
    output:
        r0 = protected("%s/{yid}/{sid}.fq.gz" % config['barn']['fqdir']),
        r1 = protected("%s/{yid}/{sid}_1.fq.gz" % config['barn']['fqdir']),
        r2 = protected("%s/{yid}/{sid}_2.fq.gz" % config['barn']['fqdir']),
    params:
        N = "{yid}.%s.{sid}" % config['fq_copy']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['fq_copy']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['fq_copy']['id']),
        j = lambda w: get_resource(w, config, 'fq_copy'),
        mem = lambda w: get_resource(w, config, 'fq_copy')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fq_copy')['ppn']
    script: "../scripts/fq_copy.py"

def fastqc_inputs(w):
    yid, sid = w.yid, w.sid
    pre = "%s/%s/%s" % (config['barn']['fqdir'], yid, sid)
    r0, r1, r2 = ["%s%s.fq.gz" % (pre, suf) for suf in ['', '_1', '_2']]
    inputs = []
    if config['y'][yid]['t'][sid]['paired']:
        inputs += [r1, r2]
    else:
        inputs += [r0]
    return inputs

rule fastqc:
    input: fastqc_inputs
    output:
        z0 = "{yid}/%s/{sid}_fastqc.zip" % config['barn']['od15'],
        z1 = "{yid}/%s/{sid}_1_fastqc.zip" % config['barn']['od15'],
        z2 = "{yid}/%s/{sid}_2_fastqc.zip" % config['barn']['od15'],
    params:
        odir = "{yid}/%s" % config['barn']['od15'],
        pre = "{yid}/%s/{sid}" % config['barn']['od15'],
        N = "{yid}.%s.{sid}" % config['fastqc']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['fastqc']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['fastqc']['id']),
        j = lambda w: get_resource(w, config, 'fastqc'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fastqc')['ppn']
    conda: "../envs/multiqc.yml"
    script: "../scripts/fastqc.py"

def multiqc_inputs(w):
    yid = w.yid
    inputs = []
    for sid in config['y'][yid]['SampleID']:
        pre = "%s/%s/%s" % (yid, config['barn']['od15'], sid)
        r0, r1, r2 = ["%s%s_fastqc.zip" % (pre, suf) for suf in ['', '_1', '_2']]
        inputs += [r0, r1, r2]
    return inputs

rule multiqc:
    input: multiqc_inputs
    output: "{yid}/%s/multiqc_data/multiqc_fastqc.txt" % config['barn']['od15']
    params:
        odir = "{yid}/%s" % config['barn']['od15'],
        N = "{yid}.%s" % config['multiqc']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['multiqc']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['multiqc']['id']),
        j = lambda w: get_resource(w, config, 'multiqc'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'multiqc')['ppn']
    conda: "../envs/multiqc.yml"
    shell: "multiqc -f -o {params.odir} {params.odir}"

def update_readlist_inputs(w):
    yid = w.yid
    idir = config['barn']['idir_sra']
    if config['y'][yid]['source'] == 'local':
        idir = config['barn']['idir_local']

    sl = "%s/%s/%s.tsv" % (config['dirh'], idir, yid),
    fqc = "%s/%s/multiqc_data/multiqc_fastqc.txt" % (yid, config['barn']['od15'])
    return dict(sl=sl, fqc=fqc)

rule update_readlist:
    input: unpack(update_readlist_inputs)
    output: "%s/%s/{yid}.tsv" % (config['dirh'], config['barn']['odir'])
    params:
        N = "{yid}.%s" % config['update_readlist']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['update_readlist']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['update_readlist']['id']),
        j = lambda w: get_resource(w, config, 'update_readlist'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'update_readlist')['ppn']
    conda: "../envs/r.yml"
    shell: "samplelist_addstat.R {input.sl} {input.fqc} {output}"

