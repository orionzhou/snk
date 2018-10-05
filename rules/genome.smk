rule fasta:
    input:
        "{genome}/08_seq_map/raw.fna"
    output:
        fna = "{genome}/10_genome.fna",
        fai = "{genome}/10_genome.fna.fai",
        chrom_bed = "{genome}/15_intervals/01.chrom.bed",
        chrom_size = "{genome}/15_intervals/01.chrom.sizes",
        gap_bed = "{genome}/15_intervals/11.gap.bed",
    params:
        wdir = lambda wildcards: "%s" % wildcards.genome,
        odir = "08_seq_map",
        gap = lambda wildcards: config['genome'][wildcards.genome]['gap'],
        prefix = lambda wildcards: config['genome'][wildcards.genome]['prefix'],
    shell:
        """
        rm -rf {output.fna}* {output.fai}*
        rm -rf {output.chrom_bed} {output.chrom_size} {output.gap_bed}
        
        mkdir -p {params.wdir}/{params.odir}
        cd {params.wdir}/{params.odir}
        rm -rf raw.fna.* renamed* map* raw.sizes
        fasta size raw.fna > raw.sizes 

        fasta rename raw.fna raw.sizes renamed.fna mapf.bed mapb.bed \
                --merge_short --gap {params.gap} --prefix_chr {params.prefix}

        fasta size renamed.fna > renamed.sizes

        chain frombed mapf.bed raw.sizes renamed.sizes > mapf.chain

        chainSwap mapf.chain mapb.chain
        
        cd ..
        ln -sf {params.odir}/renamed.fna 10_genome.fna
        cd ..
        
        samtools faidx {output.fna}
        fasta size --bed {output.fna} > {output.chrom_bed}
        cut -f1,3 {output.chrom_bed} > {output.chrom_size}
        fasta gaps {output.fna} > {output.gap_bed}
        """

rule blat:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/db.2bit" % config['blat']['odir'],
        "{genome}/21_dbs/%s/db.2bit.tile11.ooc" % config['blat']['odir'],
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['blat']['odir'])
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        faToTwoBit {input} {output[0]}
        blat {output[0]} tmp.fas tmp.out -makeOoc={output[1]}
        """

rule bwa:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/db.bwt" % config['bwa']['odir']
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['bwa']['odir'])
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        bwa index -a bwtsw -p {params.odir}/db {input}
        """

rule star:
    input:
        fna = "{genome}/10_genome.fna",
        gtf = "{genome}/50_annotation/10.gtf"
    output:
        "{genome}/21_dbs/%s/SA" % config['star']['odir']
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['star']['odir'])
    threads:
        config['star']['threads']
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        STAR --runThreadN {threads} --runMode genomeGenerate \
                --genomeDir {params.odir}/ \
                --genomeFastaFiles {input.fna} --sjdbGTFfile {input.gtf}
        """

rule gatk:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/db.fasta" % config['gatk']['odir'],
        "{genome}/21_dbs/%s/db.dict" % config['gatk']['odir'],
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['gatk']['odir'])
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        cp -f {input} {output[0]}
        gatk CreateSequenceDictionary -R {output[0]}
        samtools faidx {output[0]}
        """


