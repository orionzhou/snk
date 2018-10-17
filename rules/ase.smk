rule ase:
    input:
        "%s/{sid}.bam" % config['ase']['idir']
    output:
        protected("%s/{sid}.tsv" % config['ase']['odir'])
    log:
        "%s/ase/{sid}.log" % config['dirl']
    params:
        pre = "%s/{sid}" % config['ase']['odir'],
        gbed = config['ase']['gene_bed'],
        vbed = lambda wildcards: config['vbed'][wildcards.sid],
        extra = '',
        N = lambda w: "ase.%s" % (w.sid),
        ppn = config['ase']['ppn'],
        walltime = config['ase']['walltime'],
        mem = config['ase']['mem']
    threads: config['ase']['ppn']
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

rule merge_ase:
    input:
        expand(["%s/{sid}.tsv" % config['merge_ase']['idir']], sid = config['SampleID'])
    output:
        protected("%s/%s" % (config['dird'], config['merge_ase']['out']))
    shell:
        "merge.ase.R -o {output} {input}"


