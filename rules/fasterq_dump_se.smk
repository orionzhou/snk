rule fasterq_dump:
    input:
    output:
        "%s/{sid}.fq.gz" % config['fasterq_dump']['odir']
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

