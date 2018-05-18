rule fastqc:
    input:
        #"10.fastq/{sid}.fq.gz"
        "14.trim/{sid}.fq.gz"
    output:
        html="qc/fastqc/{sid}.html",
        zip="qc/fastqc/{sid}.zip"
    params:
        config["fastqc"]["extra"]
    wrapper:
        "0.23.1/bio/fastqc"


