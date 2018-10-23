rule fasterq_dump_pe:
    input:
    output:
        protected("%s/{sid}_1.fq.gz" % config['fasterq_dump']['odir']),
        protected("%s/{sid}_2.fq.gz" % config['fasterq_dump']['odir'])
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['fasterq_dump']['id'])
    params:
        outdir = config['fasterq_dump']['odir'],
        o1 = "%s/{sid}_1.fastq" % config['fasterq_dump']['odir'],
        o2 = "%s/{sid}_2.fastq" % config['fasterq_dump']['odir'],
        tmp = config['tmpdir'],
        N = lambda w: "%s.%s" % (config['fasterq_dump']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fasterq_dump']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fasterq_dump']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fasterq_dump')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fasterq_dump')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fasterq_dump')['mem']
    threads: config['fasterq_dump']['ppn']
    shell:
        """
        fasterq-dump --split-files -e {threads} -m {params.mem} \
                -O {params.outdir} -t {params.tmp} {wildcards.sid} \
                >{log} 2>&1

        pigz -p {threads} --fast -c {params.o1} >{output[0]}
        pigz -p {threads} --fast -c {params.o2} >{output[1]}
        rm {params.o1}
        rm {params.o2}
        """

rule fasterq_dump_se:
    input:
    output:
        protected("%s/{sid}.fq.gz" % config['fasterq_dump']['odir'])
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['fasterq_dump']['id'])
    params:
        outdir = config['fasterq_dump']['odir'],
        o1 = "%s/{sid}.fastq" % config['fasterq_dump']['odir'],
        tmp = config['tmpdir'],
        N = lambda w: "%s.%s" % (config['fasterq_dump']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fasterq_dump']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fasterq_dump']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fasterq_dump')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fasterq_dump')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fasterq_dump')['mem']
    threads: config['fasterq_dump']['ppn']
    shell:
        """
        fasterq-dump --split-files -e {threads} -m {params.mem} \
                -O {params.outdir} -t {params.tmp} {wildcards.sid} \
                >{log} 2>&1

        pigz -p {threads} --fast -c {params.o1} >{output}
        rm {params.o1}
        """

