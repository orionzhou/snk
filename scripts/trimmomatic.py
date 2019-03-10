from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config
if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    shell("""
    trimmomatic PE -threads {threads} \
    {input.r1} {input.r2} {output.r1} {output.r1u} {output.r2} {output.r2u} \
    {params.trimmer} >{log} 2>&1
    touch {output.r0}
    """)
else:
    shell("""
    trimmomatic SE -threads {threads} \
    {input} {output} \
    {params.trimmer} >{log} 2>&1
    touch {output.r1} {output.r2} {output.r1u} {output.r2u}
    """)
