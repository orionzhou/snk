rule deinterleave:
    input:
        "%s/{sid}.fq.gz" % config['deinterleave']['idir']
    output:
        r1 = "%s/{sid}_1.fq.gz" % config['deinterleave']['odir'],
        r2 = "%s/{sid}_2.fq.gz" % config['deinterleave']['odir']
    params:
        config["deinterleave"]["extra"]
    run:
        shell("zcat {input} | "
        "deinterleave_fastq.sh {output.r1} {output.r2} compress"
        )


