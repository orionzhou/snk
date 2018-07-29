rule fastq_dump:
    input:
    output:
        "%s/{sid}.fq.gz" % config['fastq_dump']['odir']
    params:
        outdir = config['fastq_dump']['odir'],
        o1 = "%s/{sid}_1.fastq.gz" % config['fastq_dump']['odir']
    run:
        shell("fastq-dump --gzip --split-files -outdir {params.outdir} {wildcards.sid}")
        shell("mv {params.o1} {output}")

