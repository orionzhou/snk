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
        "0.23.1/bio/multiqc"


