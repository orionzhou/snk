def star_extra(wildcards):
    extras = [config["star"]["extra"]]
    extras.append("--outSAMattrRGline ID:%s SM:%s" % (wildcards.sid, wildcards.sid))
    #if 'vcf' in config and wildcards.gt in config['vcf']:
    if 1 == 2:
        extras.append("--varVCFfile %s" % config['vcf'][wildcards.gt]) 
        extras.append("--waspOutputMode SAMtag")
        extras.append("--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch vG vA vW")
    else:
        extras.append("--outSAMattributes All")
    return " ".join(extras)

rule star_se:
    input: 
        "%s/%s.fq.gz" % (config['star']['idir'], sid)
    output:
        temp("%s/{sid}/Aligned.out.bam" % config['star']['odir'][0])
    log:
        "%s/star/{sid}.log" % config['dirl']
    params:
        index = config["star"]["index"],
        extra = star_extra
    threads:
        config["star"]["threads"]
    wrapper:
        "0.27.0/bio/star/align"

rule sambamba_sort:
    input:
        "%s/{sid}/Aligned.out.bam" % config['star']['odir'][0]
    output: 
        protected("%s/{sid}.bam" % config['star']['odir'][1])
    params:
        config['sambamba']['sort']['extra']
    threads:
        config['sambamba']['threads']
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input}"
 
