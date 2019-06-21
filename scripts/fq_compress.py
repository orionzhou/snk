import os.path as op
from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    if input.r1.endswith(".gz"):
        try:
            shell("ln -f %s %s" % (op.abspath(input.r1), output.r1))
        except:
            shell("cp -fL %s %s" % (op.abspath(input.r1), output.r1))
        try:
            shell("ln -f %s %s" % (op.abspath(input.r2), output.r2))
        except:
            shell("cp -fL %s %s" % (op.abspath(input.r2), output.r2))
    else:
        shell("pigz -p {threads} -c {input.r1} > {output.r1}")
        shell("pigz -p {threads} -c {input.r2} > {output.r2}")

    shell("touch {output.r0}")
else:
    if input.r0.endswith(".gz"):
        try:
            shell("ln -f %s %s" % (op.abspath(input.r0), output.r0))
        except:
            shell("cp -fL %s %s" % (op.abspath(input.r0), output.r0))
    else:
        shell("pigz -p {threads} -c {input.r0} > {output.r0}")

    shell("touch {output.r1} {output.r2}")
