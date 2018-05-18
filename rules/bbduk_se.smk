rule bbduk:
    input:
        "10.fastq/{sid}.fq.gz"
    output:
        "14.trim/{sid}.fq.gz"
    log:
        "logs/bbduk/{sid}.log"
    params:
        extra = "ref=%s %s" % 
            (config['bbduk']['refs'], config['bbduk']['extra'])
    run:
        shell("bbduk.sh "
        "in={input} "
        "out={output} "
        "{params.extra} "
        "> {log}")
