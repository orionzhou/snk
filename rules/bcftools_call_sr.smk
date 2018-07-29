rule bcftools_call:
    input:
        "%s/{sid}.bam" % config['bcftools']['idir']
    output:
        temp("%s/{sid}/{region}.bcf" % config['bcftools']['odir'])
    params:
        ref = config['bcftools']['ref'],
        extra = config['bcftools']['call']['extra']
    log:
        "%s/bcftools/{sid}/{region}.log" % config['dirl']
    threads:
        config["bcftools"]['call']["threads"]
    shell:
        """
        (bcftools mpileup -f {params.ref} -g 1 -r {wildcards.region} {input} -O u | \
        bcftools call -m -g 1 {params.extra} -O b -o {output} -) \
        2>{log}
        """

rule bcftools_concat:
    input:
        lambda wildcards: ["%s/%s/%s.bcf" % 
            (config['bcftools']['odir'], wildcards.sid, x) 
            for x in config['bcftools']['concat']['regions'].split(" ")]
    output:
        temp("%s/{sid}.bcf" % config['bcftools']['odir']),
        temp("%s/{sid}.bcf.csi" % config['bcftools']['odir'])
    params:
        input_str = lambda wildcards, input: ["%s" % x for x in input],
        extra = config['bcftools']['concat']['extra']
    log:
        "%s/bcftools/{sid}.log" % config['dirl']
    threads:
        config["bcftools"]['concat']["threads"]
    shell:
        "bcftools concat {params.extra} -O b -o {output[0]} {input} >{log} 2>&1"
        "bcftools index {output[0]}"

rule bcftools_merge:
    input:
        expand(["%s/{sid}.bcf" % config['bcftools']['odir']], sid = config['t']['SampleID'])
    output:
        protected(config['bcftools']['outfile'])
    params:
        ref = config['bcftools']['ref'],
        gvcf = "%s/all.bcf" % config['bcftools']['odir'],
        extra = config['bcftools']['merge']['extra']
    log:
        "%s/bcftools.log" % config['dirl']
    threads:
        config["bcftools"]['merge']["threads"]
    shell:
        """
        (bcftools merge --threads {threads} -g {params.ref} {input} -O u | \
        bcftools view --threads {threads} -m 2 -O v -o {output} -) \
        2>{log}
        """
