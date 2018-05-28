def star_inputs(wildcards):
    sid = wildcards.sid
    inputs = dict()
    if config['paired']:
        inputs['fq1'] = "14.trim/%s_1.PE.fq.gz" % sid
        inputs['fq2'] = "14.trim/%s_2.PE.fq.gz" % sid
    else:
        inputs['fq1'] = "14.trim/%s.fq.gz" % sid
    return inputs

def star_extra(wildcards):
    extras = [config["star"]["extra"]]
    extras.append("--outSAMattrRGline ID:%s SM:%s" % (wildcards.sid, wildcards.sid))
    if 'vcf' in config:
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
        "21.bam.raw/{sid}_{gt}/Aligned.out.bam"
    log:
        "logs/star/{sid}_{gt}.log"
    params:
        index = config["star"]["index"],
        outprefix = "21.bam.raw/{sid}_{gt}/",
        extra = star_extra
    threads:
        config["star"]["threads"]
    wrapper:
        "0.23.1/bio/star/align"

