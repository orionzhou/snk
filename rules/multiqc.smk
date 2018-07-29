rule multiqc:
    input:
        multiqc_inputs
    output:
        "%s/multiqc.html" % config['dirq']
    params:
        outdir = lambda wildcards, output: op.dirname(output[0]),
        outfile = lambda wildcards, output: op.basename(output[0]),
        extra = config["multiqc"]["extra"]
    log:
        "%s/multiqc.log" % config['dirl']
    shell:
        "multiqc {params.extra} --force "
        "-o {params.outdir} "
        "-n {params.outfile} "
        "{input} "
        ">{log} 2>&1"

