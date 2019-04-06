rule fastqc:
    input:
        "%s/{sid}.fq.gz" % config['fastqc']['idir']
    output:
        html="%s/fastqc/{sid}_fastqc.html" % config['dirq'],
        zip="%s/fastqc/{sid}_fastqc.zip" % config['dirq']
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['fastqc']['id'])
    params:
        extra = '',
        N = lambda w: "%s.%s" % (config['ase']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirj'], config['fastqc']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirj'], config['fastqc']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fastqc')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fastqc')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fastqc')['mem']
    conda: "../envs/job.yml"
    wrapper:
        "0.27.0/bio/fastqc"


