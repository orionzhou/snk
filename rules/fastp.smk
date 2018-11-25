rule fastp_se:
    input:
        fq = "%s/{sid}.fq.gz" % config['fastp']['idir']
    output:
        fq = "%s/{sid}.fq.gz" % config['fastp']['odir'],
        json = "%s/{sid}.se.json" % config['fastp']['odir'],
        html = "%s/{sid}.se.html" % config['fastp']['odir'],
    params:
        N = lambda w: "%s.%s" % (config['fastp']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fastp']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fastp']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fastp')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fastp')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fastp')['mem']
    threads: config['fastp']['ppn']
    shell:
        """
        fastp --thread {threads} \
        -i {input.fq} -o {output.fq} \
        -j {output.json} -h {output.html}
        """

rule fastp_pe:
    input:
        r1 ="%s/{sid}_1.fq.gz" % config['fastp']['idir'],
        r2 ="%s/{sid}_2.fq.gz" % config['fastp']['idir']
    output:
        r1 = "%s/{sid}_1.fq.gz" % config['fastp']['odir'],
        r2 = "%s/{sid}_2.fq.gz" % config['fastp']['odir'],
        json = "%s/{sid}.pe.json" % config['fastp']['odir'],
        html = "%s/{sid}.pe.html" % config['fastp']['odir'],
    params:
        N = lambda w: "%s.%s" % (config['fastp']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fastp']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fastp']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fastp')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fastp')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fastp')['mem']
    threads: config['fastp']['ppn']
    shell:
        """
        fastp --thread {threads} \
        -i {input.r1} -I {input.r2} \
        -o {output.r1} -O {output.r2} \
        -j {output.json} -h {output.html}
        """
