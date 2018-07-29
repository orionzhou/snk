rule bcftools_gtcheck:
    input:
        "%s/{sid}.bcf" % config['bcftools']['odir']
    output:
        "%s/{sid}.txt" % config['bcftools']['gtcheck']['odir']
    params:
        vcf = config['bcftools']['gtcheck']['vcf'],
        extra = config['bcftools']['gtcheck']['extra']
    shell:
        """
        bcftools gtcheck {params.extra} -g {params.vcf} {input} > {output}
        """

