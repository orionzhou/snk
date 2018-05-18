def multiqc_inputs(wildcards):
    inputs = ["24.featurecounts/01.txt"]
    t = config['t']
    for i in range(len(t)):
        sid, gt = t['sid'][i], t['genotype'][i]
        if config['paired']:
            inputs.append("qc/fastqc/%s_1.PE.zip" % sid)
            inputs.append("qc/fastqc/%s_1.SE.zip" % sid)
            inputs.append("qc/fastqc/%s_2.PE.zip" % sid)
            inputs.append("qc/fastqc/%s_2.SE.zip" % sid)
        else:
            inputs.append("qc/fastqc/%s.zip" % sid)
        inputs.append("21.bam.raw/%s_%s/Aligned.out.bam" % (sid, gt))
    return inputs

rule multiqc:
    input:
        multiqc_inputs
    output:
        "qc/multiqc.html"
    params:
        config["multiqc"]["extra"]
    log:
        "logs/multiqc.log"
    wrapper:
        "0.23.0/bio/multiqc"


