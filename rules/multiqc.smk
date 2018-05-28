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
        "0.23.1/bio/multiqc"


