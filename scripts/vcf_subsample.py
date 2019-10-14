import pandas as pd
from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config


def infer_subsample_fraction(fs, nsites=1000000):
    sl = pd.read_csv(fs, sep="\t", header=0)
    total_sites = 0
    for i in range(len(sl)):
        if sl['key'][i] == 'number of records':
            total_sites = int(sl['value'][i])
            break
    return nsites / total_sites

fraction = infer_subsample_fraction(input['stat'], params['nsites']),

shell("""
    gatk --java-options "-Xmx{params.mem}G" SelectVariants \
        --tmp-dir {params.tmp} \
        --use-jdk-deflater --use-jdk-inflater \
        -R {params.refn} \
        -V {input.vcf} -O {output} -fraction {fraction}
    bcftools stats -s - {output} > {params.stat}
    """)

