def hisat2_inputs(w):
    sid = w.sid
    idir = config['hisat2']['idir']
    inputs = dict()
    if config['t'][sid]['paired']:
        inputs['fq1'] = "%s/%s_1.fq.gz" % (idir, sid)
        inputs['fq2'] = "%s/%s_2.fq.gz" % (idir, sid)
    else:
        inputs['fq'] = "%s/%s.fq.gz" % (idir, sid)
    return inputs

def hisat2_input_str(w):
    sid = w.sid
    idir = config['hisat2']['idir']
    input_str = ''
    if config['t'][sid]['paired']:
        fq1 = "%s/%s_1.fq.gz" % (idir, sid)
        fq2 = "%s/%s_2.fq.gz" % (idir, sid)
        input_str = "-1 %s -2 %s" % (fq1, fq2)
    else:
        fq = "%s/%s.fq.gz" % (idir, sid)
        input_str = "-U %s" % fq
    return input_str

def hisat2_extra(w):
    extras = [config["hisat2"]["extra"]]
    extras.append("--rg-id %s --rg SM:%s" % (w.sid, w.sid))
    return " ".join(extras)

rule hisat2:
    input:
        unpack(hisat2_inputs)
    output:
        temp("%s/{sid}.bam" % config['hisat2']['odir1']),
        "%s/{sid}.txt" % config['hisat2']['odir1']
    params:
        index = config[config['reference']]["hisat2"],
        input_str = hisat2_input_str,
        extra = hisat2_extra,
        N = lambda w: "%s.%s" % (config['hisat2']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['hisat2']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['hisat2']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'hisat2')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'hisat2')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'hisat2')['mem']
    threads: config['hisat2']['ppn']
    shell:
        """
        hisat2 {params.extra} --threads {threads} \
            -x {params.index} {params.input_str} \
            --summary-file {output[1]} \
            | samtools view -Sbh -o {output[0]} -
        """

rule sambamba_sort:
    input:
        "%s/{sid}.bam" % config['hisat2']['odir1']
    output:
        "%s/{sid}.bam" % config['hisat2']['odir2']
    params:
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra']),
        N = lambda w: "%s.%s" % (config['sambamba']['sort']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['sambamba']['sort']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['sambamba']['sort']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['mem']
    threads: config['sambamba']['ppn']
    shell:
        "sambamba sort {params.extra} -t {threads} -o {output} {input}"

