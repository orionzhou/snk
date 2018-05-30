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
    if 'vcf' in config and wildcards.gt in config['vcf']:
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
        "%s/{sid}_{gt}/Aligned.out.bam" % config['star']['odir'][0]
    log:
        "%s/star/{sid}_{gt}.log" % config['dirl']
    params:
        index = config["star"]["index"],
        outprefix = "%s/{sid}_{gt}/" % config['star']['odir'][0],
        extra = star_extra
    threads:
        config["star"]["threads"]
    wrapper:
        "0.23.1/bio/star/align"

rule sambamba_sort:
    input:
        "%s/{pre}/Aligned.out.bam" % config['star']['odir'][0]
    output: 
        "%s/{pre}.bam" % config['star']['odir'][1]
    params:
        config['sambamba']['sort']['extra']
    threads:
        config['sambamba']['threads']
    wrapper:
        "0.23.1/bio/sambamba/sort"
 
