from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config
if params.mapper == 'bwa':
    shell("""
    sambamba view -S -f bam -t {threads} {input} -o {output[0]}
    sambamba sort {params.extra} -t {threads} -o {output[0]} {params.tmp_bam}
    rm {params.tmp_bam}
    """)
else:
    shell("sambamba sort {params.extra} -t {threads} -o {output[0]} {input}")
