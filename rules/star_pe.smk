def star_inputs(wildcards):
    sid = wildcards.sid
    inputs = dict()
    if config['paired']:
        inputs['fq1'] = "%s/%s_1.fq.gz" % (config['star']['idir'], sid)
        inputs['fq2'] = "%s/%s_2.fq.gz" % (config['star']['idir'], sid)
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

rule star_pe:
    input: 
        unpack(star_inputs)
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

rule star_pe_unpaired:
    input:
        fq1 = [
            "%s/{sid}_1.unpaired.fq.gz" % (config['star']['idir']),
            "%s/{sid}_2.unpaired.fq.gz" % (config['star']['idir'])
        ]
    output:
        temp("%s/{sid}_unpaired/Aligned.out.bam" % config['star']['odir'][0])
    log:
        "%s/star/{sid}.unpaired.log" % config['dirl']
    params:
        index = config["star"]["index"],
        extra = star_extra
    threads:
        config["star"]["threads"]
    wrapper:
        "0.27.0/bio/star/align"

rule sambamba_merge:
    input:
        "%s/{sid}/Aligned.out.bam" % config['star']['odir'][0],
        "%s/{sid}_unpaired/Aligned.out.bam" % config['star']['odir'][0],
    output: 
        temp("%s/{sid}.unsorted.bam" % config['star']['odir'][1])
    params:
        config['sambamba']['merge']['extra']
    threads:
        config['sambamba']['threads']
    shell:
        "sambamba merge {params} -t {threads} -o {output} {input}"
 
rule sambamba_sort:
    input:
        "%s/{sid}.unsorted.bam" % config['star']['odir'][1]
    output: 
        protected("%s/{sid}.bam" % config['star']['odir'][1])
    params:
        config['sambamba']['sort']['extra']
    threads:
        config['sambamba']['threads']
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input}"
