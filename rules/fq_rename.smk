rule copy_fq:
    input:
        "%s/{sid}.fq.gz" % config['copy_fq']['idir']
    output:
        protected("%s/{sid}.fq.gz" % config['copy_fq']['odir'])
    shell:
        "cp -f {input} {output}"

rule gzip_fq:
    input:
        "%s/{sid}.fq" % config['copy_fq']['idir']
    output:
        protected("%s/{sid}.fq.gz" % config['copy_fq']['odir'])
    shell:
        "gzip -c {input} > {output}"

