rule trimmomatic:
    input:
        "10.fastq/{sid}.fq.gz"
    output:
        "14.trim/{sid}.fq.gz"
    log:
        "logs/trimmomatic/{sid}.log"
    params:
        extra = "-threads %s ILLUMINACLIP:%s:2:30:10:8:no "
        "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 " % 
        (config["trimmomatic"]["threads"], 
        config['trimmomatic']['adapter_se'])
    wrapper:
        "0.23.1/bio/trimmomatic/se"
