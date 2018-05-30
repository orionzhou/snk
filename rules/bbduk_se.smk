rule bbduk:
    input:
        "%s/{sid}.fq.gz" % config['bbduk']['idir']
    output:
        "%s/{sid}.fq.gz" % config['bbduk']['odir']
    log:
        "%s/bbduk/{sid}.log" % config['dirl']
    params:
        cmd = config['bbduk']['cmd'],
        extra = "ref=%s %s" % 
            (','.join(config['bbduk']['refs']), config['bbduk']['extra']),
    run:
        shell("{params.cmd} "
        "in={input} "
        "out={output} "
        "{params.extra} "
        "> {log}")
