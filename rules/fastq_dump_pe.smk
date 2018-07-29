rule fastq_dump:
    input:
    output:
        "%s/{sid}_1.fq.gz" % config['fastq_dump']['odir'],
        "%s/{sid}_2.fq.gz" % config['fastq_dump']['odir']
    params:
        outdir = config['fastq_dump']['odir'],
        o1 = "%s/{sid}_1.fastq.gz" % config['fastq_dump']['odir'],
        o2 = "%s/{sid}_2.fastq.gz" % config['fastq_dump']['odir']
    run:
        shell("fastq-dump --gzip --split-files -outdir {params.outdir} {wildcards.sid}")
        shell("mv {params.o1} {output[0]}")
        shell("mv {params.o2} {output[1]}")

