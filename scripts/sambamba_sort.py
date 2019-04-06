from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config
if params.mapper == 'bwa':
    shell("""
    sambamba view -S -f bam -t {threads} {input} -o {params.bam0}
    sambamba sort {params.extra} -t {threads} -o {output[0]} {params.bam0}
    rm {params.bam0}
    """)
else:
    shell("sambamba sort {params.extra} -t {threads} -o {output[0]} {input}")
