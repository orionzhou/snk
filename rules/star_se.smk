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

rule star:
    input: 
        "%s/{sid}.fq.gz" % (config['star']['idir'])
    output:
        temp("%s/{sid}/Aligned.out.bam" % config['star']['odir'][0]),
        protected("%s/{sid}/Log.final.out" % config['star']['odir'][0])
    log:
        "%s/star/{sid}.log" % config['dirl']
    params:
        index = config["star"]["index"],
        input_str = lambda wildcards, input: " ".join(input),
        outprefix = "%s/{sid}/" % config['star']['odir'][0],
        readcmd = lambda wildcards, input: "--readFilesCommand zcat" if input[0].endswith(".gz") else "",
        extra = star_extra
    threads:
        config["star"]["threads"]
    shell:
        """
        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix} \
        --outStd Log \
        >{log} 2>&1
        """

rule sambamba_sort:
    input:
        "%s/{sid}/Aligned.out.bam" % config['star']['odir'][0]
    output: 
        protected("%s/{sid}.bam" % config['star']['odir'][1])
    params:
        "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra'])
    threads:
        config['sambamba']['threads']
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input}"
 
