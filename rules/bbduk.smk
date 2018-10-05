rule bbduk_se:
    input:
        "%s/{sid}.fq.gz" % config['bbduk']['idir']
    output:
        "%s/{sid}.fq.gz" % config['bbduk']['odir'],
        "%s/bbduk_se/{sid}.log" % config['dirl']
    params:
        cmd = config['bbduk']['cmd'],
        extra = "ref=%s %s" % 
            (','.join(config['bbduk']['refs']), config['bbduk']['extra']),
    shell:
        "{params.cmd} in={input} out={output[0]} {params.extra} stats={output[1]}"
