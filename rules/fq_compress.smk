rule fq_compress_se:
    input:
        i1 = "%s/{sid}.fq" % config['fq_compress']['idir']
    output:
        o1 = protected("%s/{sid}.fq.gz" % config['fq_compress']['odir'])
    threads: config['fq_compress']['threads']
    shell:
        "pigz -p {threads} -c {input.i1} > {output.o1}"

rule fq_compress_pe:
    input:
        i1 = "%s/{sid}_1.fq" % config['fq_compress']['idir'],
        i2 = "%s/{sid}_2.fq" % config['fq_compress']['idir']
    output:
        o1 = protected("%s/{sid}_1.fq.gz" % config['fq_compress']['odir']),
        o2 = protected("%s/{sid}_2.fq.gz" % config['fq_compress']['odir'])
    threads: config['fq_compress']['threads']
    shell:
        """
        pigz -p {threads} -c {input.i1} > {output.o1}
        pigz -p {threads} -c {input.i2} > {output.o2}
        """

