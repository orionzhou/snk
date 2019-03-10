from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config
if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    shell("""
    fasterq-dump --split-files -e {threads} -m {params.mem} \
            -O {params.odir} -t {params.tmp} {wildcards.sid}

    pigz -p {threads} --fast -c {params.o1} >{output.r1}
    pigz -p {threads} --fast -c {params.o2} >{output.r2}
    rm {params.o1} {params.o2}
    touch {output.r0}
    """)
else:
    shell("""
    fasterq-dump --split-files -e {threads} -m {params.mem} \
            -O {params.odir} -t {params.tmp} {wildcards.sid}

    pigz -p {threads} --fast -c {params.o0} >{output.r0}
    rm {params.o0}
    touch {output.r1} {output.r2}
    """)
