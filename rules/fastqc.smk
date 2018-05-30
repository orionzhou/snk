rule fastqc:
    input:
        "%s/{sid}.fq.gz" % config['fastqc']['idir']
    output:
        html="%s/fastqc/{sid}_fastqc.html" % config['dirq'],
        zip="%s/fastqc/{sid}_fastqc.zip" % config['dirq']
    params:
        config["fastqc"]["extra"]
    wrapper:
        "0.23.1/bio/fastqc"


