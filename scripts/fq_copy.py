import os.path as op
from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

try:
    shell("ln -f %s %s" % (op.abspath(input.r0), output.r0))
except:
    shell("cp -fL %s %s" % (op.abspath(input.r0), output.r0))
try:
    shell("ln -f %s %s" % (op.abspath(input.r1), output.r1))
except:
    shell("cp -fL %s %s" % (op.abspath(input.r1), output.r1))
try:
    shell("ln -f %s %s" % (op.abspath(input.r2), output.r2))
except:
    shell("cp -fL %s %s" % (op.abspath(input.r2), output.r2))
