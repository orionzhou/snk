rule trimmomatic:
    input:
        "10.fastq/{sid}.fq.gz"
    output:
        "14.trim/{sid}.fq.gz"
    log:
        "logs/trimmomatic/{sid}.log"
    params:
        extra = "-threads %s" % config["trimmomatic"]["threads"],
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_se'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"]
    wrapper:
        "0.23.1/bio/trimmomatic/se"
