rule ase_old:
    input: "%s/{sid}.bam" % config['ase']['idir']
    output: protected("%s/{sid}.tsv" % config['ase']['odir'])
    params:
        pre = "%s/{sid}" % config['ase']['odir'],
        gbed = config['ase']['gene_bed'],
        vbed = lambda w: config['vbed'][w.sid],
        extra = '',
        N = lambda w: "%s.%s" % (config['ase']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirj'], config['ase']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirj'], config['ase']['id'], w.sid),
        j = lambda w: get_resource(w, config, 'ase'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase')['ppn']
    shell:
        """
        mkdir -p {params.pre}
        ase.py bam2bed {input} {params.pre}.1.bed
        sort -T {params.pre} -k1,1 -k2,2n {params.pre}.1.bed -o {params.pre}.2.sorted.bed
        intersectBed -wa -wb -a {params.pre}.2.sorted.bed -b {params.vbed} > {params.pre}.3.bed
        sort -T {params.pre} -k4,4 -k1,1 -k2,2n {params.pre}.3.bed > {params.pre}.4.sorted.bed
        ase.py bed_prep {params.pre}.4.sorted.bed {params.pre}.5.tsv {params.pre}.6.bed
        sort -T {params.pre} -k1,1 -k2,2n {params.pre}.6.bed -o {params.pre}.7.sorted.bed
        intersectBed -wa -wb -a {params.gbed} -b {params.pre}.7.sorted.bed > {params.pre}.8.bed
        ase.py bed_summarise {params.pre}.5.tsv {params.pre}.8.bed {params.pre}.tsv
        rm {params.pre}.[1-8].*
        rm -rf {params.pre}
        """


rule merge_ase:
    input:
        expand(["{yid}/%s/{sid}.tsv" % config['merge_ase']['idir']], sid = config['SampleID'])
    output: protected("%s/%s" % (config['oid'], config['merge_ase']['out']))
    shell:
        "merge.ase.R -o {output} {input}"


