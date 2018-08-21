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

rule star_se:
    input: 
        "%s/{sid}.fq.gz" % config['star']['idir']
    output:
        temp("%s/{sid}/Aligned.out.bam" % config['star']['odir1']),
        protected("%s/{sid}/Log.final.out" % config['star']['odir1'])
    log:
        "%s/star/{sid}.log" % config['dirl']
    params:
        index = config["star"]["index"],
        input_str = lambda wildcards, input: " ".join(input),
        outprefix = "%s/{sid}/" % config['star']['odir1'],
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

rule star_pe:
    input: 
        fq1 = "%s/{sid}_1.fq.gz" % config['star']['idir'],
        fq2 = "%s/{sid}_2.fq.gz" % config['star']['idir'],
        fq1u = "%s/{sid}_1.unpaired.fq.gz" % config['star']['idir'],
        fq2u = "%s/{sid}_2.unpaired.fq.gz" % config['star']['idir']
    output:
        temp("%s/{sid}_p/Aligned.out.bam" % config['star']['odir1']),
        temp("%s/{sid}_u/Aligned.out.bam" % config['star']['odir1']),
        protected("%s/{sid}_p/Log.final.out" % config['star']['odir1']),
        protected("%s/{sid}_u/Log.final.out" % config['star']['odir1'])
    log:
        "%s/star/{sid}.log" % config['dirl']
    params:
        index = config["star"]["index"],
        input_str_p = lambda wildcards, input: "%s %s" % (input.fq1, input.fq2),
        input_str_u = lambda wildcards, input: "%s,%s" % (input.fq1u, input.fq2u),
        outprefix_p = "%s/{sid}_p/" % config['star']['odir1'],
        outprefix_u = "%s/{sid}_u/" % config['star']['odir1'],
        readcmd = lambda wildcards, input: "--readFilesCommand zcat" if input.fq1.endswith(".gz") else "",
        extra = star_extra
    threads:
        config["star"]["threads"]
    shell:
        """
        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str_p} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix_p} \
        --outStd Log \
        >{log} 2>&1

        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str_u} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix_u} \
        --outStd Log \
        >{log} 2>&1
        """

def sambamba_sort_inputs(wildcards):
    sid = wildcards.sid
    odir1 = config['star']['odir1']
    inputs = dict()
    if config['t'][sid]['paired']:
        inputs['bam_p'] = "%s/%s_p/Aligned.out.bam" % (odir1, sid)
        inputs['bam_u'] = "%s/%s_u/Aligned.out.bam" % (odir1, sid)
    else:
        inputs['bam'] = "%s/%s/Aligned.out.bam" % (odir1, sid)
    return inputs

rule sambamba_sort:
    input:
        unpack(sambamba_sort_inputs)
    output: 
        protected("%s/{sid}.bam" % config['star']['odir2'])
    params:
        sorted_p = "%s/{sid}_p/Aligned.out.sorted.bam" % (config['star']['odir1']),
        sorted_u = "%s/{sid}_u/Aligned.out.sorted.bam" % (config['star']['odir1']),
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra'])
    threads:
        config['sambamba']['threads']
    run:
        if config['t'][wildcards.sid]['paired']:
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_p} {input.bam_p}")
            shell("sambamba sort {params.extra} -t {threads} -o {params.sorted_u} {input.bam_u}")
            shell("sambamba merge -t {threads} {output} {params.sorted_p} {params.sorted_u}")
            shell("rm {params.sorted_p} {params.sorted_u}")
        else:
            shell("sambamba sort {params.extra} -t {threads} -o {output} {input.bam}")

