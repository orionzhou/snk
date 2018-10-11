rule deinterleave:
    input:
        "%s/{sid}.fq.gz" % config['deinterleave']['idir']
    output:
        r1 = protected("%s/{sid}_1.fq.gz" % config['deinterleave']['odir']),
        r2 = protected("%s/{sid}_2.fq.gz" % config['deinterleave']['odir'])
    params:
        config["deinterleave"]["extra"]
    threads:
        config['deinterleave']['threads']
    shell:
        "zcat {input} | "
        "deinterleave_fastq.sh {output.r1} {output.r2} compress"

