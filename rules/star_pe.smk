def star_inputs(wildcards):
    sid = wildcards.sid
    inputs = dict()
    if config['paired']:
        inputs['fq1'] = "%s/%s_1.fq.gz" % (config['star']['idir'], sid)
        inputs['fq2'] = "%s/%s_2.fq.gz" % (config['star']['idir'], sid)
        inputs['fq1u'] = "%s/%s_1.unpaired.fq.gz" % (config['star']['idir'], sid)
        inputs['fq2u'] = "%s/%s_2.unpaired.fq.gz" % (config['star']['idir'], sid)
    else:
        inputs['fq1'] = "%s/%s.fq.gz" % (config['star']['idir'], sid)
    return inputs

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
        unpack(star_inputs)
    output:
        temp("%s/{sid}/Aligned.out.bam" % config['star']['odir'][0]),
        temp("%s/{sid}_unpaired/Aligned.out.bam" % config['star']['odir'][0]),
        protected("%s/{sid}/Log.final.out" % config['star']['odir'][0]),
        protected("%s/{sid}_unpaired/Log.final.out" % config['star']['odir'][0])
    log:
        "%s/star/{sid}.log" % config['dirl']
    params:
        index = config["star"]["index"],
        input_str_paired = lambda wildcards, input: "%s %s" % (input.fq1, input.fq2),
        input_str_unpaired = lambda wildcards, input: "%s,%s" % (input.fq1u, input.fq2u),
        outprefix_paired = "%s/{sid}/" % config['star']['odir'][0],
        outprefix_unpaired = "%s/{sid}_unpaired/" % config['star']['odir'][0],
        readcmd = lambda wildcards, input: "--readFilesCommand zcat" if input.fq1.endswith(".gz") else "",
        extra = star_extra
    threads:
        config["star"]["threads"]
    shell:
        """
        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str_paired} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix_paired} \
        --outStd Log \
        >{log} 2>&1

        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str_unpaired} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix_unpaired} \
        --outStd Log \
        >{log} 2>&1
        """

rule sambamba_sort:
    input:
        "%s/{sid}/Aligned.out.bam" % config['star']['odir'][0],
        "%s/{sid}_unpaired/Aligned.out.bam" % config['star']['odir'][0]
    output: 
        temp("%s/{sid}/Aligned.out.sorted.bam" % config['star']['odir'][0]),
        temp("%s/{sid}_unpaired/Aligned.out.sorted.bam" % config['star']['odir'][0])
    params:
        "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra'])
    threads:
        config['sambamba']['threads']
    shell:
        """
        sambamba sort {params} -t {threads} -o {output[0]} {input[0]}
        sambamba sort {params} -t {threads} -o {output[1]} {input[1]}
        """

rule sambamba_merge:
    input:
        "%s/{sid}/Aligned.out.sorted.bam" % config['star']['odir'][0],
        "%s/{sid}_unpaired/Aligned.out.sorted.bam" % config['star']['odir'][0]
    output: 
        protected("%s/{sid}.bam" % config['star']['odir'][1])
    params:
        config['sambamba']['merge']['extra']
    threads:
        config['sambamba']['threads']
    shell:
        "sambamba merge {params} -t {threads} {output} {input}"
 
