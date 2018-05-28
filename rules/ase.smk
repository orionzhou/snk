def ase_gbed(wildcards):
    gbed = config['annotation']['bed']
    if 'gene_bed' in config['ase']:
        gbed = config['ase']['gene_bed']
    return gbed

rule ase:
    input:
        "22.bam/{sid}_{gt}.bam"
    output:
        "26.ase/{sid}_{gt}.tsv"
    params:
        pre = "26.ase/{sid}_{gt}.bam",
        gbed = ase_gbed,
        vbed = lambda wildcards: config['vbed'][wildcards.gt],
        extra = config['ase']['extra'],
    run:
        shell("mkdir -p {params.pre}"),
        shell("bam2bed.py {input} {params.pre}.1.bed"),
        shell("sort -T {params.pre} -k1,1 -k2,2n {params.pre}.1.bed -o {params.pre}.2.sorted.bed"),
        shell("intersectBed -wa -wb -a {params.pre}.2.sorted.bed -b {params.vbed} > {params.pre}.3.bed"),
        shell("sort -T {params.pre} -k4,4 -k1,1 -k2,2n {params.pre}.3.bed > {params.pre}.4.sorted.bed"),
        shell("bed.ase.py {params.pre}.4.sorted.bed {params.pre}.5.tsv {params.pre}.6.bed"),
        shell("sort -T {params.pre} -k1,1 -k2,2n {params.pre}.6.bed -o {params.pre}.7.sorted.bed"),
        shell("intersectBed -wa -wb -a {params.gbed} -b {params.pre}.7.sorted.bed > {params.pre}.8.bed"),
        shell("bed.ase.sum.py {params.pre}.5.tsv {params.pre}.8.bed {params.pre}.tsv"),
        shell("rm {params.pre}.[1-8].*"),
        shell("rm -rf {params.pre}")

