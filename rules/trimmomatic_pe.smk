rule trimmomatic:
    input:
        r1 ="%s/{sid}_1.fq.gz" % config['trimmomatic']['idir'],
        r2 ="%s/{sid}_2.fq.gz" % config['trimmomatic']['idir']
    output:
        r1 = "%s/{sid}_1.fq.gz" % config['trimmomatic']['odir'],
        r2 = "%s/{sid}_2.fq.gz" % config['trimmomatic']['odir'],
        r1_unpaired = "%s/{sid}_1.unpaired.fq.gz" % config['trimmomatic']['odir'],
        r2_unpaired = "%s/{sid}_2.unpaired.fq.gz" % config['trimmomatic']['odir']
    log:
        "%s/trimmomatic/{sid}.log" % config['dirl']
    params:
        extra = "-threads %s" % config["trimmomatic"]["threads"],
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_pe'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"]
    wrapper:
        "0.23.1/bio/trimmomatic/pe"

