import os
from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

shell("mkdir -p {params.odir}")
os.chdir(params.odir)
for genome in config['ortho_genomes']:
    fi = config['g'][genome]['annotation']['lfaa']
    fo = "%s.faa" % genome
    shell("ln -sf {fi} {fo}")

shell("orthofinder -f . -t {threads} -p . -a 4")
shell("rm -rf output")
shell("mv OrthoFinder/Results_* output")
shell("rm -rf OrthoFinder")
shell("rm *.faa")
