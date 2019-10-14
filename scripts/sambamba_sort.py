import os.path as op
from snakemake import shell
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

sort_tag = ''
if params.get('sort', '') == 'byname':
    sort_tag = '-n'

fos = []
for fi in input:
    fname = op.basename(fi)
    if params.get('mapper','') == 'star':
        fname = op.basename(op.dirname(fi))
    fo = op.join(params.get('odir'), fname)
    shell("sambamba sort {sort_tag} --tmpdir={params.tmp_dir} -t {threads} -o {fo} {fi}")
    fos.append(fo)

if len(fos) == 1:
    shell("mv {fos[0]} {output[0]}")
else:
    bam_str = ' '.join(fos)
    shell("sambamba merge -t {threads} {output[0]} {bam_str}")
    bam_str = ' '.join([x + '*' for x in fos])
    shell("rm {bam_str}")

if params.get('sort', '') != 'byname':
    shell("samtools index {output[0]}")
