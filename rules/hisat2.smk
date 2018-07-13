def hisat2_inputs(wildcards):
    sid = wildcards.sid
    inputs = []
    if config['paired']:
        inputs.append("%s/%s_1.fq.gz" % (config['hisat2']['idir'], sid))
        inputs.append("%s/%s_2.fq.gz" % (config['hisat2']['idir'], sid))
    else:
        inputs.append("%s/%s.fq.gz" % (config['hisat2']['idir'], sid))
    return { 'reads': inputs }

def hisat2_extra(wildcards):
    extras = [config["hisat2"]["extra"]]
    extras.append("--rg-id %s --rg SM:%s" % (wildcards.sid, wildcards.sid))
    return " ".join(extras)

rule hisat2:
    input: 
        unpack(hisat2_inputs)
    output:
        "%s/{sid}_{gt}.bam" % config['hisat2']['odir'][0]
    log:
        "%s/hisat2/{sid}_{gt}.log" % config['dirl']
    params:
        idx = config["hisat2"]["index"],
        extra = hisat2_extra
    threads:
        config["hisat2"]["threads"]
    wrapper:
        "0.27.0/bio/hisat2"

rule sambamba_sort:
    input:
        "%s/{pre}.bam" % config['hisat2']['odir'][0]
    output: 
        "%s/{pre}.bam" % config['hisat2']['odir'][1]
    params:
        config['sambamba']['sort']['extra']
    threads:
        config['sambamba']['threads']
    wrapper:
        "0.27.0/bio/sambamba/sort"
 
