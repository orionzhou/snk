rule deinterleave:
    input:
        "%s/{sid}.fq.gz" % config['deinterleave']['idir']
    output:
        r1 = protected("%s/{sid}_1.fq.gz" % config['deinterleave']['odir']),
        r2 = protected("%s/{sid}_2.fq.gz" % config['deinterleave']['odir'])
    params:
        extra = config["deinterleave"]["extra"],
        N = lambda w: "%s.%s" % (config['deinterleave']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['deinterleave']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['deinterleave']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'deinterleave')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'deinterleave')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'deinterleave')['mem']
    threads: config['deinterleave']['ppn']
    shell:
        "zcat {input} | "
        "deinterleave_fastq.sh {output.r1} {output.r2} {threads} compress"

