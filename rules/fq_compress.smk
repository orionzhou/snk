rule fq_compress_se:
    input:
        i1 = "%s/{sid}.fq" % config['fq_compress']['idir']
    output:
        o1 = protected("%s/{sid}.fq.gz" % config['fq_compress']['odir'])
    params:
        N = lambda w: "%s.%s" % (config['fq_compress']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fq_compress']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fq_compress']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['mem']
    threads: config['fq_compress']['ppn']
    shell:
        "pigz -p {threads} -c {input.i1} > {output.o1}"

rule fq_compress_pe:
    input:
        i1 = "%s/{sid}_1.fq" % config['fq_compress']['idir'],
        i2 = "%s/{sid}_2.fq" % config['fq_compress']['idir']
    output:
        o1 = protected("%s/{sid}_1.fq.gz" % config['fq_compress']['odir']),
        o2 = protected("%s/{sid}_2.fq.gz" % config['fq_compress']['odir'])
    params:
        N = lambda w: "%s.%s" % (config['fq_compress']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fq_compress']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fq_compress']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fq_compress')['mem']
    threads: config['fq_compress']['ppn']
    shell:
        """
        pigz -p {threads} -c {input.i1} > {output.o1}
        pigz -p {threads} -c {input.i2} > {output.o2}
        """

