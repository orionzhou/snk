rule fasta:
    input: "{genome}/download/raw.fna"
    output:
        fna = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        fai = "{genome}/%s/%s.fai" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        chrom_size = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['chrom_size']),
        chrom_bed = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['chrom_bed']),
        gap = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['gap']),
        fchain = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['fchain']),
        bchain = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['bchain']),
    params:
        wdir = lambda w: "%s" % w.genome,
        odir = "08_seq_map",
        gap = lambda w: config['x'][w.genome]['gap'],
        prefix = lambda w: config['x'][w.genome]['prefix'],
        N = "{genome}.%s" % config['fasta']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['fasta']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['fasta']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fasta')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fasta')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fasta')['mem']
    threads: config['fasta']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {output.fna}* {output.fai}*
        rm -rf {output.chrom_bed} {output.chrom_size} {output.gap}

        mkdir -p {params.wdir}/{params.odir}
        cd {params.wdir}/{params.odir}
        rm -rf raw.fna.* renamed* map* raw.sizes
        ln -sf ../download/raw.fna raw.fna
        fasta.py size raw.fna > raw.sizes

        fasta.py rename raw.fna raw.sizes renamed.fna mapf.bed mapb.bed \
                --merge_short --gap {params.gap} --prefix_chr {params.prefix}

        fasta.py size renamed.fna > renamed.sizes

        chain.py fromBed mapf.bed raw.sizes renamed.sizes > mapf.chain

        chainSwap mapf.chain mapb.chain

        cd ..
        ln -sf {params.odir}/renamed.fna 10_genome.fna
        cd ..

        samtools faidx {output.fna}
        fasta.py size --bed {output.fna} > {output.chrom_bed}
        cut -f1,3 {output.chrom_bed} > {output.chrom_size}
        fasta.py gaps {output.fna} > {output.gap}
        """

rule blat_index:
    input:
        "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref'])
    output:
        "{genome}/%s/%s" % (config['db']['blat']['xdir'], config['db']['blat']['x.2bit']),
        "{genome}/%s/%s" % (config['db']['blat']['xdir'], config['db']['blat']['x.ooc']),
    params:
        odir = "{genome}/%s" % config['db']['blat']['xdir'],
        N = "{genome}.%s" % config['blat_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['blat_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['blat_index']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'blat_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'blat_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'blat_index')['mem']
    threads: config['blat_index']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        faToTwoBit {input} {output[0]}
        blat {output[0]} tmp.fas tmp.out -makeOoc={output[1]}
        """

rule blast_index:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
    output:
        "{genome}/%s/%s" % (config['db']['blastn']['xdir'], config['db']['blastn']['xout']),
        "{genome}/%s/%s" % (config['db']['blastp']['xdir'], config['db']['blastp']['xout']),
    params:
        odir1 = "{genome}/%s" % config['db']['blastn']['xdir'],
        odir2 = "{genome}/%s" % config['db']['blastp']['xdir'],
        N = "{genome}.%s" % config['blast_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['blast_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['blast_index']['id']),
    threads: config['blast_index']['ppn']
    conda: "../envs/blast.yml"
    shell:
        """
        makeblastdb -dbtype nucl -in {input.fna} -title db -out {params.odir1}/db
        makeblastdb -dbtype prot -in {input.faa} -title db -out {params.odir2}/db
        """

rule bwa_index:
    input:
        "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref'])
    output:
        "{genome}/%s/%s" % (config['db']['bwa']['xdir'], config['db']['bwa']['xout'])
    params:
        odir = "{genome}/%s" % config['db']['bwa']['xdir'],
        N = "{genome}.%s" % config['bwa_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['bwa_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['bwa_index']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bwa_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bwa_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bwa_index')['mem']
    threads: config['bwa_index']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        bwa index -a bwtsw -p {params.odir}/db {input}
        """

rule bismark_index:
    input:
        "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref'])
    output:
        "{genome}/%s/%s" % (config['db']['bismark']['xdir'], config['db']['bismark']['xout'])
    params:
        odir = "{genome}/%s" % config['db']['bismark']['xdir'],
        parallel = lambda w, resources: int(resources.ppn / 2),
        N = "{genome}.%s" % config['bismark_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['bismark_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['bismark_index']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bismark_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bismark_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bismark_index')['mem']
    threads: config['bismark_index']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        cd {params.odir}
        ln -sf ../../10_genome.fna db.fa
        bismark_genome_preparation --bowtie2 .
        """

rule star_index:
    input:
        fna = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        gtf = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gtf'])
    output:
        "{genome}/%s/%s" % (config['db']['star']['xdir'], config['db']['star']['xout'])
    params:
        odir = "{genome}/%s" % config['db']['star']['xdir'],
        N = "{genome}.%s" % config['star_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['star_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['star_index']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'star_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'star_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'star_index')['mem']
    threads: config['star_index']['ppn']
    conda: "../envs/work.yml"
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
        "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref'])
    output:
        "{genome}/%s/%s" % (config['db']['gatk']['xdir'], config['db']['gatk']['xref']),
        "{genome}/%s/%s" % (config['db']['gatk']['xdir'], config['db']['gatk']['xref.dict']),
    params:
        odir = "{genome}/%s" % config['db']['gatk']['xdir'],
        N = "{genome}.%s" % config['gatk_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['gatk_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['gatk_index']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk_index')['mem']
    threads: config['gatk_index']['ppn']
    conda: "../envs/work.yml"
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
        fna = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        gtf = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gtf'])
    output:
        "{genome}/%s/%s" % (config['db']['hisat2']['xdir'], config['db']['hisat2']['xout'])
    params:
        odir = "{genome}/%s" % config['db']['hisat2']['xdir'],
        N = "{genome}.%s" % config['hisat2_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['hisat2_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['hisat2_index']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'hisat2_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'hisat2_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'hisat2_index')['mem']
    threads: config['hisat2_index']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        hisat2_extract_exons.py {input.gtf} > {params.odir}/db.exon
        hisat2_extract_splice_sites.py {input.gtf} > {params.odir}/db.ss
        hisat2-build -p {threads} --ss {params.odir}/db.ss \
                --exon {params.odir}/db.exon {input.fna} {params.odir}/db
        """

rule snpeff_index:
    input:
        fna = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        gff = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lgff']),
    output:
        "{genome}/%s/%s" % (config['db']['snpeff']['xdir'], config['db']['snpeff']['xcfg']),
        "{genome}/%s/{genome}/%s" % (config['db']['snpeff']['xdir'], config['db']['snpeff']['xout'])
    params:
        fna = lambda w, input: op.abspath(input.fna),
        gff = lambda w, input: op.abspath(input.gff),
        odir = "{genome}/%s" % config['db']['snpeff']['xdir'],
        odir1 = "{genome}/%s/{genome}" % config['db']['snpeff']['xdir'],
        N = "{genome}.%s" % config['snpeff_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['snpeff_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['snpeff_index']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'snpeff_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'snpeff_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'snpeff_index')['mem']
    threads: config['snpeff_index']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        mkdir -p {params.odir1}
        echo 'data.dir = .' > {output[0]}
        echo '{wildcards.genome}.genome : Zea mays' >> {output[0]}
        ln -sf {params.fna} {params.odir1}/sequences.fa
        ln -sf {params.gff} {params.odir1}/genes.gff
        snpEff -Xmx{params.mem}G build -c {output[0]} -gff3 -v {wildcards.genome}
        """

