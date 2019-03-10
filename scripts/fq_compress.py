from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config
if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    if input.r1.endswith(".gz"):
        shell("""
        ln -sf {input.r1} {output.r1}
        ln -sf {input.r2} {output.r2}
        touch {output.r0}
        """)
    else:
        shell("""
        pigz -p {threads} -c {input.r1} > {output.r1}
        pigz -p {threads} -c {input.r2} > {output.r2}
        touch {output.r0}
        """)
else:
    if input.r0.endswith(".gz"):
        shell("""
        ln -sf {input.r0} {output.r0}
        touch {output.r1} {output.r2}
        """)
    else:
        shell("""
        pigz -p {threads} -c {input.r0} > {output.r0}
        touch {output.r1} {output.r2}
        """)
