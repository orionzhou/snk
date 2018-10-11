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
        gap = lambda wildcards: config['genomes'][wildcards.genome]['gap'],
        prefix = lambda wildcards: config['genomes'][wildcards.genome]['prefix'],
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

rule blat_index:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/db.2bit" % config['blat']['xdir'],
        "{genome}/21_dbs/%s/db.2bit.tile11.ooc" % config['blat']['xdir'],
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['blat']['xdir'])
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        faToTwoBit {input} {output[0]}
        blat {output[0]} tmp.fas tmp.out -makeOoc={output[1]}
        """

rule bwa_index:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/db.bwt" % config['bwa']['xdir']
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['bwa']['xdir'])
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        bwa index -a bwtsw -p {params.odir}/db {input}
        """

rule star_index:
    input:
        fna = "{genome}/10_genome.fna",
        gtf = "{genome}/50_annotation/10.gtf"
    output:
        "{genome}/21_dbs/%s/SA" % config['star']['xdir']
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['star']['xdir'])
    threads:
        config['star']['xthreads']
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        STAR --runThreadN {threads} --runMode genomeGenerate \
                --genomeDir {params.odir}/ \
                --genomeFastaFiles {input.fna} --sjdbGTFfile {input.gtf}
        """

rule gatk_index:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/db.fasta" % config['gatk']['xdir'],
        "{genome}/21_dbs/%s/db.dict" % config['gatk']['xdir'],
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['gatk']['xdir'])
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        cp -f {input} {output[0]}
        gatk CreateSequenceDictionary -R {output[0]}
        samtools faidx {output[0]}
        """

rule hisat2_index:
    input:
        fna = "{genome}/10_genome.fna",
        gtf = "{genome}/50_annotation/10.gtf"
    output:
        "{genome}/21_dbs/%s/db.1.ht2" % config['hisat2']['xdir']
    params:
        odir = lambda wildcards: "%s/21_dbs/%s" % (wildcards.genome, config['hisat2']['xdir'])
    threads:
        config['hisat2']['xthreads']
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        hisat2_extract_exons.py {input.gtf} > {params.odir}/db.exon
        hisat2_extract_splice_sites.py {input.gtf} > {params.odir}/db.ss
        hisat2-build -p {threads} --ss {params.odir}/db.ss \
                --exon {params.odir}/db.exon {input.fna} {params.odir}/db
        """


