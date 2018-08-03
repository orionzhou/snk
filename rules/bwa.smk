rule bwa_se:
    input: 
        "%s/{sid}.fq.gz" % config['bwa']['idir']
    output:
        temp("%s/{sid}.bam" % config['bwa']['odir'])
    params:
        index = config["bwa"]["index"],
        extra = config['bwa']['extra']
    threads:
        config["bwa"]["threads"]
    shell:
        """
        (bwa mem -t {threads} {params.index} {input} | \
        samtools view -Sbh -o {output} -) \
        >{log} 2>&1
        """

rule bwa_pe:
    input: 
        fq1 = "%s/{sid}_1.fq.gz" % config['bwa']['idir'],
        fq2 = "%s/{sid}_2.fq.gz" % config['bwa']['idir'],
        fq1u = "%s/{sid}_1.unpaired.fq.gz" % config['bwa']['idir'],
        fq2u = "%s/{sid}_2.unpaired.fq.gz" % config['bwa']['idir']
    output:
        temp("%s/{sid}_p.bam" % config['bwa']['odir']),
        temp("%s/{sid}_u1.bam" % config['bwa']['odir']),
        temp("%s/{sid}_u2.bam" % config['bwa']['odir'])
    log:
        "%s/bwa/{sid}.log" % config['dirl']
    params:
        index = config["bwa"]["index"],
        extra = config['bwa']['extra']
    threads:
        config["bwa"]["threads"]
    shell:
        """
        (bwa mem -t {threads} {params.index} {input.fq1} {input.fq2} | \
        samtools view -Sbh -o {output[0]} -) \
        >{log} 2>&1
        (bwa mem -t {threads} {params.index} {input.fq1u} | \
        samtools view -Sbh -o {output[1]} -) \
        >>{log} 2>&1
        (bwa mem -t {threads} {params.index} {input.fq2u} | \
        samtools view -Sbh -o {output[2]} -) \
        >>{log} 2>&1
        """

def sambamba_sort_inputs(wildcards):
    sid = wildcards.sid
    odir = config['bwa']['odir']
    inputs = dict()
    if config['t'][sid]['paired']:
        inputs['bam_p'] = "%s/%s_p.bam" % (odir, sid)
        inputs['bam_u1'] = "%s/%s_u1.bam" % (odir, sid)
        inputs['bam_u2'] = "%s/%s_u2.bam" % (odir, sid)
    else:
        inputs['bam'] = "%s/%s.bam" % (odir, sid)
    return inputs

rule sambamba_sort:
    input:
        unpack(sambamba_sort_inputs)
    output: 
        protected("%s/{sid}.sorted.bam" % config['bwa']['odir'])
    params:
        sorted_p = "%s/{sid}_p.bam" % (config['bwa']['odir']),
        sorted_u1 = "%s/{sid}_u1.bam" % (config['bwa']['odir']),
        sorted_u2 = "%s/{sid}_u2.bam" % (config['bwa']['odir']),
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra'])
    threads:
        config['sambamba']['threads']
    run:
        if config['t'][wildcards.sid]['paired']:
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_p} {input.bam_p}")
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_u1} {input.bam_u1}")
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_u2} {input.bam_u2}")
            shell("sambamba merge -t {threads} {output} {params.sorted_p} {params.sorted_u1} {params.sorted_u2}")
        else:
            shell("sambamba sort {params.extra} -t {threads} -o {output} {input.bam}")

rule sambamba_flagstat:
    input:
        "%s/{sid}.sorted.bam" % config['bwa']['odir']
    output:
        protected("%s/{sid}.txt" % config['bwa']['odir'])
    params:
        extra = config['sambamba']['flagstat']['extra']
    threads:
        config['sambamba']['threads']
    shell:
        "sambamba flagstat {params.extra} -t {threads} {input} > {output}"
