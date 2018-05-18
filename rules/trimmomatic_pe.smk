rule trimmomatic:
    input:
        r1 ="10.fastq/{sid}_1.fq.gz",
        r2 ="10.fastq/{sid}_2.fq.gz"
    output:
        r1p = "14.trim/{sid}_1.PE.fq.gz",
        r2p = "14.trim/{sid}_1.SE.fq.gz",
        r1u = "14.trim/{sid}_2.PE.fq.gz",
        r2u = "14.trim/{sid}_2.SE.fq.gz"
    log:
        "logs/trimmomatic/{sid}.log"
    params:
        extra = "-threads %s ILLUMINACLIP:%s:2:30:10:8:no "
        "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 " % 
        (config["trimmomatic"]["threads"], 
        config['trimmomatic']['adapter_pe'])
    wrapper:
        "0.23.1/bio/trimmomatic/pe"

