rule bcftools_call:
    input:
        bams = expand(["%s/{sid}.bam" % config['callvnt']['idir']], 
                sid = config['SampleID']),
        bais = expand(["%s/{sid}.bam.bai" % config['callvnt']['idir']], 
                sid = config['SampleID'])
    output:
        protected("%s/{rid}.vcf.gz" % config['callvnt']['odir2']),
        protected("%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2'])
    params:
        ref = config['bcftools']['ref'],
        region = lambda wildcards: config['regions'][wildcards.rid],
        extra = config['bcftools']['call']['extra']
    threads:
        config["bcftools"]['call']["threads"]
    shell:
        """
        bcftools mpileup -f {params.ref} -r {params.region} \
                {input.bams} -O u | \
                bcftools call -c {params.extra} -O z -o {output[0]} -
        
        bcftools index -t {output[0]}
        """

rule bcftools_concat:
    input:
        vcfs = expand(["%s/{rid}.vcf.gz" % config['callvnt']['odir2']], 
                rid = list(config['regions'].keys())),
        tbis = expand(["%s/{rid}.vcf.gz.tbi" % config['callvnt']['odir2']], 
                rid = list(config['regions'].keys()))
    output:
        vcf = protected("%s" % config['callvnt']['outfile']),
        tbi = protected("%s.tbi" % config['callvnt']['outfile'])
    params:
        extra = config['bcftools']['concat']['extra']
    threads:
        config["bcftools"]['concat']["threads"]
    shell:
        """
        bcftools concat {params.extra} -O z -o {output.vcf} \
                {input.vcfs}
        bcftools index -t {output.vcf}
        """

