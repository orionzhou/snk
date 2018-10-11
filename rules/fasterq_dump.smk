rule fasterq_dump_pe:
    input:
    output:
        protected("%s/{sid}_1.fq.gz" % config['fasterq_dump']['odir']),
        protected("%s/{sid}_2.fq.gz" % config['fasterq_dump']['odir'])
    params:
        outdir = config['fasterq_dump']['odir'],
        o1 = "%s/{sid}_1.fastq" % config['fasterq_dump']['odir'],
        o2 = "%s/{sid}_2.fastq" % config['fasterq_dump']['odir'],
        mem = config['fasterq_dump']['mem'],
        tmp = config['tmpdir']
    threads:
        config["fasterq_dump"]["threads"]
    log:
        "%s/fasterq_dump/{sid}.log" % config['dirl']
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
    params:
        outdir = config['fasterq_dump']['odir'],
        o1 = "%s/{sid}.fastq" % config['fasterq_dump']['odir'],
        mem = config['fasterq_dump']['mem'],
        tmp = config['tmpdir']
    threads:
        config["fasterq_dump"]["threads"]
    log:
        "%s/fasterq_dump/{sid}.log" % config['dirl']
    shell:
        """
        fasterq-dump --split-files -e {threads} -m {params.mem} \
                -O {params.outdir} -t {params.tmp} {wildcards.sid} \
                >{log} 2>&1

        pigz -p {threads} --fast -c {params.o1} >{output}
        rm {params.o1}
        """
