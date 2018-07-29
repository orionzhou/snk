rule bcftools_call:
    input:
        "%s/{sid}.bam" % config['bcftools']['idir']
    output:
        protected("%s/{sid}.bcf" % config['bcftools']['odir']),
        protected("%s/{sid}.bcf.csi" % config['bcftools']['odir'])
    params:
        ref = config['bcftools']['ref'],
        extra = config['bcftools']['call']['extra']
    log:
        "%s/bcftools/{sid}.log" % config['dirl']
    threads:
        config["bcftools"]['call']["threads"]
    shell:
        """
        (bcftools mpileup -f {params.ref} -g 1 {input} -O u | \
        bcftools call -m -g 1 {params.extra} -O b -o {output[0]} -) \
        2>{log}
        
        bcftools index {output[0]}
        """

rule bcftools_merge:
    input:
        bcfs = expand(["%s/{sid}.bcf" % config['bcftools']['odir']], sid = config['t']['SampleID']),
        csis = expand(["%s/{sid}.bcf.csi" % config['bcftools']['odir']], sid = config['t']['SampleID'])
    output:
        temp("%s/{region}.bcf" % config['bcftools']['odir']),
        temp("%s/{region}.bcf.csi" % config['bcftools']['odir'])
    params:
        ref = config['bcftools']['ref'],
        reg = lambda wildcards: config['regions'][wildcards.region],
        extra = config['bcftools']['merge']['extra']
    log:
        "%s/bcftools/{region}.log" % config['dirl']
    threads:
        config["bcftools"]['merge']["threads"]
    shell:
        """
        (bcftools merge --threads {threads} -g {params.ref} -r {params.reg} {input.bcfs} -O u | \
        bcftools view --threads {threads} -m 2 -O b -o {output[0]} -) \
        2>{log}

        bcftools index {output[0]}
        """

rule bcftools_concat:
    input:
        bcfs = expand(["%s/{region}.bcf" % config['bcftools']['odir']], 
                       region = list(config['regions'].keys())),
        csis = expand(["%s/{region}.bcf.csi" % config['bcftools']['odir']],
                       region = list(config['regions'].keys()))
    output:
        protected(config['bcftools']['outfile'])
    params:
        extra = config['bcftools']['concat']['extra']
    log:
        "%s/bcftools.log" % config['dirl']
    threads:
        config["bcftools"]['concat']["threads"]
    shell:
        "bcftools concat {params.extra} -O b -o {output} {input.bcfs} >{log} 2>&1"
        "bcftools index {output}"

