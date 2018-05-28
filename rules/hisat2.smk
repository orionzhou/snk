def hisat2_inputs(wildcards):
    sid = wildcards.sid
    inputs = []
    if config['paired']:
        inputs.append("14.trim/%s_1.PE.fq.gz" % sid)
        inputs.append("14.trim/%s_2.PE.fq.gz" % sid)
    else:
        inputs.append("14.trim/%s.fq.gz" % sid)
    return { 'reads': inputs }

def hisat2_extra(wildcards):
    extras = [config["hisat2"]["extra"]]
    extras.append("--rg-id %s --rg SM:%s" % (wildcards.sid, wildcards.sid))
    return " ".join(extras)

rule hisat2:
    input: 
        unpack(hisat2_inputs)
    output:
        "21.bam.raw/{sid}_{gt}.bam"
    log:
        "logs/hisat2/{sid}_{gt}.log"
    params:
        idx = config["hisat2"]["index"],
        extra = hisat2_extra
    threads:
        config["hisat2"]["threads"]
    wrapper:
        "0.23.1/bio/hisat2"

