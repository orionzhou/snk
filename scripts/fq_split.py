import os.path as op
from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

yid, sid = wildcards.yid, wildcards.sid
npart = config['y'][wildcards.yid]['t'][wildcards.sid]['npart']
parts = config['y'][wildcards.yid]['t'][wildcards.sid]['parts']

fho = open(output[0], 'w')
for part in config['y'][wildcards.yid]['t'][wildcards.sid]['parts']:
    print(part, file=fho)
fho.close()

line_per_file = config['trimming']['part_size'] * 4

pre = "%s/{sid}" % params.odir
pre0 = "%s_" % pre
pre1 = "%s_1_" % pre
pre2 = "%s_2_" % pre
if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
    shell("zcat %s | split -a 1 -l %d - %s" % (input.r1, line_per_file, pre1))
    shell("zcat %s | split -a 1 -l %d - %s" % (input.r2, line_per_file, pre2))
else:
    shell("zcat %s | split -a 1 -l %d - %s" % (input.r0, line_per_file, pre0))

for part in parts:
    i_r0 = "%s%s" % (pre0, part)
    i_r1 = "%s%s" % (pre1, part)
    i_r2 = "%s%s" % (pre2, part)
    o_r0 = "%s_%s.fq" % (pre, part)
    o_r1 = "%s_%s_1.fq" % (pre, part)
    o_r2 = "%s_%s_2.fq" % (pre, part)
    if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
        shell("mv %s %s" % (i_r1, o_r1))
        shell("mv %s %s" % (i_r2, o_r2))
    else:
        shell("mv %s %s" % (i_r0, o_r0))
        # shell("pigz -p {threads} -c %s > %s" % (i_r1, o_r1))
        # shell("pigz -p {threads} -c %s > %s" % (i_r2, o_r2))
        # shell("pigz -p {threads} -c %s > %s" % (i_r0, o_r0))
