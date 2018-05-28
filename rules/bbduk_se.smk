rule bbduk:
    input:
        "10.fastq/{sid}.fq.gz"
    output:
        "14.trim/{sid}.fq.gz"
    log:
        "logs/bbduk/{sid}.log"
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
