rule trimmomatic_se:
    input:
        "%s/{sid}.fq.gz" % config['trimmomatic']['idir']
    output:
        "%s/{sid}.fq.gz" % config['trimmomatic']['odir']
    log:
        "%s/%s_se/{sid}.log" % (config['dirl'], config['trimmomatic']['id'])
    params:
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_se'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"],
        N = lambda w: "%s.%s" % (config['trimmomatic']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['trimmomatic']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['trimmomatic']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['mem']
    threads: config['trimmomatic']['ppn']
    shell:
        "trimmomatic SE -threads {threads} "
        "{input} {output} "
        "{params.trimmer} "
        ">{log} 2>&1"

rule trimmomatic_pe:
    input:
        r1 ="%s/{sid}_1.fq.gz" % config['trimmomatic']['idir'],
        r2 ="%s/{sid}_2.fq.gz" % config['trimmomatic']['idir']
    output:
        r1 = "%s/{sid}_1.fq.gz" % config['trimmomatic']['odir'],
        r2 = "%s/{sid}_2.fq.gz" % config['trimmomatic']['odir'],
        r1_unpaired = "%s/{sid}_1.unpaired.fq.gz" % config['trimmomatic']['odir'],
        r2_unpaired = "%s/{sid}_2.unpaired.fq.gz" % config['trimmomatic']['odir']
    log:
        "%s/%s_pe/{sid}.log" % (config['dirl'], config['trimmomatic']['id'])
    params:
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_pe'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"],
        N = lambda w: "%s.%s" % (config['trimmomatic']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['trimmomatic']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['trimmomatic']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'trimmomatic')['mem']
    threads: config['trimmomatic']['ppn']
    shell:
        "trimmomatic PE -threads {threads} "
        "{input.r1} {input.r2} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} "
        ">{log} 2>&1"

