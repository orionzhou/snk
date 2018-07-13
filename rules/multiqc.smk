rule multiqc:
    input:
        multiqc_inputs
    output:
        "%s/multiqc.html" % config['dirq']
    params:
        config["multiqc"]["extra"]
    log:
        "%s/multiqc.log" % config['dirl']
    wrapper:
        "0.27.0/bio/multiqc"


