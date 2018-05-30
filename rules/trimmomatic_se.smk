rule trimmomatic:
    input:
        "%s/{sid}.fq.gz" % config['trimmomatic']['idir']
    output:
        "%s/{sid}.fq.gz" % config['trimmomatic']['odir']
    log:
        "%s/trimmomatic/{sid}.log" % config['dirl']
    params:
        extra = "-threads %s" % config["trimmomatic"]["threads"],
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_se'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"]
    wrapper:
        "0.23.1/bio/trimmomatic/se"
