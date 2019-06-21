from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

shell("mkdir -p {params.odir}")
shell("fastqc --threads {threads} --noextract --format fastq -o {params.odir} {input}")

if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    shell("""
    rm {params.pre}_1_fastqc.html
    rm {params.pre}_2_fastqc.html
    touch {output.z0}
    """)
else:
    shell("""
    rm {params.pre}_fastqc.html
    touch {output.z1}
    touch {output.z2}
    """)
