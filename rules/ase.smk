def ase_gbed(wildcards):
    gbed = config['annotation']['bed']
    if 'gene_bed' in config['ase']:
        gbed = config['ase']['gene_bed']
    return gbed

rule ase:
    input:
        "%s/{sid}.bam" % config['ase']['idir']
    output:
        protected("%s/{sid}.tsv" % config['ase']['odir'])
    params:
        pre = "%s/{sid}" % config['ase']['odir'],
        gbed = ase_gbed,
        vbed = lambda wildcards: config['vbed'][wildcards.gt],
        extra = config['ase']['extra'],
    log:
        "%s/ase/{sid}.log" % config['dirl']
    shell:
        """
        mkdir -p {params.pre}
        bam2bed.py {input} {params.pre}.1.bed
        sort -T {params.pre} -k1,1 -k2,2n {params.pre}.1.bed -o {params.pre}.2.sorted.bed
        intersectBed -wa -wb -a {params.pre}.2.sorted.bed -b {params.vbed} > {params.pre}.3.bed
        sort -T {params.pre} -k4,4 -k1,1 -k2,2n {params.pre}.3.bed > {params.pre}.4.sorted.bed
        bed.ase.py {params.pre}.4.sorted.bed {params.pre}.5.tsv {params.pre}.6.bed
        sort -T {params.pre} -k1,1 -k2,2n {params.pre}.6.bed -o {params.pre}.7.sorted.bed
        intersectBed -wa -wb -a {params.gbed} -b {params.pre}.7.sorted.bed > {params.pre}.8.bed
        bed.ase.sum.py {params.pre}.5.tsv {params.pre}.8.bed {params.pre}.tsv
        rm {params.pre}.[1-8].*
        rm -rf {params.pre}
        """

