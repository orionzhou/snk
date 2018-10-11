rule deinterleave:
    input:
        "%s/{sid}.fq.gz" % config['deinterleave']['idir']
    output:
        r1 = protected("%s/{sid}_1.fq.gz" % config['deinterleave']['odir']),
        r2 = protected("%s/{sid}_2.fq.gz" % config['deinterleave']['odir'])
    params:
        extra = config["deinterleave"]["extra"],
        ppn = config['deinterleave']['ppn'],
        walltime = config['deinterleave']['walltime'],
        mem = config['deinterleave']['mem']
    threads: config['deinterleave']['ppn']
    shell:
        "zcat {input} | "
        "deinterleave_fastq.sh {output.r1} {output.r2} {threads} compress"

