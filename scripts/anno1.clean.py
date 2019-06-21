from snakemake import shell
input, output, params, threads, w, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

genome = w.genome
params.hybrid = config['x'][genome]['hybrid']
params.prefix = config['x'][genome]['prefix']
if params.hybrid:
    shell("cat {input.gff1} {input.gff2} > {output}")
else:
    shell("""
        gff.py fix --opt {params.fixopt} {input.gff} > {params.wdir}/01.fixed.gff
        liftOver -gff {params.wdir}/01.fixed.gff {input.chain} \
                {params.wdir}/02.lifted.gff {params.wdir}/unmapped
        ln -sf 02.lifted.gff {output}
    """)
