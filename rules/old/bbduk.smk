rule bbduk_se:
    input:
        "%s/{sid}.fq.gz" % config['bbduk']['idir']
    output:
        "%s/{sid}.fq.gz" % config['bbduk']['odir'],
        "%s/%s_se/{sid}.log" % (config['dirl'], config['bbduk']['id'])
    params:
        cmd = config['bbduk']['cmd'],
        extra = "ref=%s %s" % 
            (','.join(config['bbduk']['refs']), config['bbduk']['extra']),
        N = lambda w: "%s.%s" % (config['bbduk']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bbduk']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bbduk']['id'], w.sid),
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
