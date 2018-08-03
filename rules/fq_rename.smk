rule copy_fq_se:
    input:
        i1 = "%s/{sid}.fq.gz" % config['fq_rename']['idir']
    output:
        o1 = protected("%s/{sid}.fq.gz" % config['fq_rename']['odir'])
    shell:
        "cp -f {input.i1} {output.o1}"

rule copy_fq_pe:
    input:
        i1 = "%s/{sid}_1.fq.gz" % config['fq_rename']['idir'],
        i2 = "%s/{sid}_2.fq.gz" % config['fq_rename']['idir']
    output:
        o1 = protected("%s/{sid}_1.fq.gz" % config['fq_rename']['odir']),
        o2 = protected("%s/{sid}_2.fq.gz" % config['fq_rename']['odir'])
    shell:
        """
        cp -f {input.i1} {output.o1}
        cp -f {input.i2} {output.o2}
        """

rule gzip_fq_se:
    input:
        i1 = "%s/{sid}.fq" % config['fq_rename']['idir']
    output:
        o1 = protected("%s/{sid}.fq.gz" % config['fq_rename']['odir'])
    shell:
        "gzip -c {input.i1} > {output.o1}"

rule gzip_fq_pe:
    input:
        i1 = "%s/{sid}.1.fq" % config['fq_rename']['idir'],
        i2 = "%s/{sid}.2.fq" % config['fq_rename']['idir']
    output:
        o1 = protected("%s/{sid}.1.fq.gz" % config['fq_rename']['odir']),
        o2 = protected("%s/{sid}.2.fq.gz" % config['fq_rename']['odir'])
    shell:
        """
        gzip -c {input.i1} > {output.o1}
        gzip -c {input.i2} > {output.o2}
        """

