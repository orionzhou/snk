def bwa_extra(w):
    extras = [config["bwa"]["extra"]]
    sm = w.sid
    if 'Genotype' in config['t'][w.sid]:
        sm = config['t'][w.sid]['Genotype']
    pl = 'ILLUMINA'
    if config['readtype'] == 'solid':
        pl = 'SOLID'
    extras.append("-R '@RG\\tID:%s\\tSM:%s\\tPL:%s'" % (w.sid, sm, pl))
    return " ".join(extras)

rule bwa_se:
    input: 
        "%s/{sid}.fq.gz" % config['bwa']['idir']
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['bwa']['id'])
    output:
        temp("%s/{sid}.sam" % config['bwa']['odir1'])
    params:
        index = config[config['reference']]["bwa"],
        extra = bwa_extra,
        N = lambda w: "%s.%s" % (config['bwa']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bwa']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bwa']['id'], w.sid),
        q = lambda w, resources: resources.q,
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'bwa')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'bwa')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'bwa')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'bwa')['mem']
    threads: config['bwa']['ppn']
    shell:
        """
        bwa mem -t {threads} {params.index} {params.extra} {input} \
                >{output} 2>>{log}
        """

rule bwa_pe:
    input: 
        fq1 = "%s/{sid}_1.fq.gz" % config['bwa']['idir'],
        fq2 = "%s/{sid}_2.fq.gz" % config['bwa']['idir'],
        fq1u = "%s/{sid}_1.unpaired.fq.gz" % config['bwa']['idir'],
        fq2u = "%s/{sid}_2.unpaired.fq.gz" % config['bwa']['idir']
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['bwa']['id'])
    output:
        temp("%s/{sid}_p.sam" % config['bwa']['odir1']),
        temp("%s/{sid}_u1.sam" % config['bwa']['odir1']),
        temp("%s/{sid}_u2.sam" % config['bwa']['odir1'])
    params:
        index = config[config['reference']]["bwa"],
        extra = bwa_extra,
        N = lambda w: "%s.%s" % (config['bwa']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bwa']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bwa']['id'], w.sid),
        q = lambda w, resources: resources.q,
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'bwa')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'bwa')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'bwa')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'bwa')['mem']
    threads: config['bwa']['ppn']
    shell:
        """
        bwa mem -t {threads} {params.index} {params.extra} {input.fq1} {input.fq2} \
        >{output[0]} 2>>{log}
        bwa mem -t {threads} {params.index} {params.extra} {input.fq1u} \
        >{output[1]} 2>>{log}
        bwa mem -t {threads} {params.index} {params.extra} {input.fq2u} \
        >{output[2]} 2>>{log}
        """

def sambamba_sort_inputs(wildcards):
    sid = wildcards.sid
    odir = config['bwa']['odir1']
    inputs = dict()
    if config['t'][sid]['paired']:
        inputs['sam_p'] = "%s/%s_p.sam" % (odir, sid)
        inputs['sam_u1'] = "%s/%s_u1.sam" % (odir, sid)
        inputs['sam_u2'] = "%s/%s_u2.sam" % (odir, sid)
    else:
        inputs['sam'] = "%s/%s.sam" % (odir, sid)
    return inputs

rule sambamba_sort:
    input:
        unpack(sambamba_sort_inputs)
    output: 
        "%s/{sid}.bam" % config['bwa']['odir2'],
        "%s/{sid}.bam.bai" % config['bwa']['odir2']
    params:
        bam = "%s/{sid}.bam" % config['bwa']['odir1'],
        bam_p = "%s/{sid}_p.bam" % config['bwa']['odir1'], 
        bam_u1 = "%s/{sid}_u1.bam" % config['bwa']['odir1'], 
        bam_u2 = "%s/{sid}_u2.bam" % config['bwa']['odir1'], 
        sorted_p = "%s/{sid}_p.sorted.bam" % config['bwa']['odir1'], 
        sorted_u1 = "%s/{sid}_u1.sorted.bam" % config['bwa']['odir1'], 
        sorted_u2 = "%s/{sid}_u2.sorted.bam" % config['bwa']['odir1'], 
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra']),
        N = lambda w: "%s.%s" % (config['sambamba']['sort']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['sambamba']['sort']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['sambamba']['sort']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'sambamba', 'sort')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'sambamba', 'sort')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'sambamba', 'sort')['mem']
    threads: config['sambamba']['ppn']
    run:
        if config['t'][wildcards.sid]['paired']:
            shell("sambamba view -S -f bam -t {threads} {input.sam_p} -o {params.bam_p}")
            shell("sambamba view -S -f bam -t {threads} {input.sam_u1} -o {params.bam_u1}")
            shell("sambamba view -S -f bam -t {threads} {input.sam_u2} -o {params.bam_u2}")
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_p} {params.bam_p}")
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_u1} {params.bam_u1}")
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_u2} {params.bam_u2}")
            shell("sambamba merge -t {threads} {output[0]} {params.sorted_p} {params.sorted_u1} {params.sorted_u2}")
            shell("rm {params.bam_p} {params.bam_u1} {params.bam_u2} {params.sorted_p}* {params.sorted_u1}* {params.sorted_u2}*")
        else:
            shell("sambamba view -S -f bam -t {threads} {input.sam} -o {params.bam}")
            shell("sambamba sort {params.extra} -t {threads} -o {output[0]} {params.bam}")
