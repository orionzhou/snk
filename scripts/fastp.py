from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config
if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    shell("""
    fastp --thread {threads} \
    -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
    -j {output.json} -h {output.html}
    touch {output.r0}
    """)
else:
    shell("""
    fastp --thread {threads} \
    -i {input.r0} -o {output.r0} \
    -j {output.json} -h {output.html}
    touch {output.r1} {output.r2}
    """)
