def fasta_input(w):
    genome = w.genome
    if config['x'][genome]['hybrid']:
        genome1, genome2 = genome.split('x')
        return ["%s/%s/%s" % (g, config['db']['fasta']['xdir'], config['db']['fasta']['ref']) for g in (genome1, genome2)]
    else:
        return "%s/download/raw.fna" % genome

rule fasta:
    input: fasta_input
    output:
        fna = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        fai = "{genome}/%s/%s.fai" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
        chrom_size = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['chrom_size']),
        chrom_bed = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['chrom_bed']),
        gap = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['gap']),
        fchain = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['fchain']),
        bchain = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['bchain']),
    params:
        wdir = "{genome}",
        odir = "08_seq_map",
        opt = lambda w: w.genome,
        N = "{genome}.%s" % config['fasta']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['fasta']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['fasta']['id']),
        j = lambda w: get_resource(w, config, 'fasta'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fasta')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/make_fasta.py"

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
        j = lambda w: get_resource(w, config, 'blat_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'blat_index')['ppn']
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
        j = lambda w: get_resource(w, config, 'blast_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'blast_index')['ppn']
    conda: "../envs/blast.yml"
    shell:
        """
        makeblastdb -dbtype nucl -in {input.fna} -title db -out {params.odir1}/db
        makeblastdb -dbtype prot -in {input.faa} -title db -out {params.odir2}/db
        """

rule last_index:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
    output:
        "{genome}/%s/%s" % (config['db']['lastn']['xdir'], config['db']['lastn']['xout']),
        "{genome}/%s/%s" % (config['db']['lastp']['xdir'], config['db']['lastp']['xout']),
    params:
        odir1 = "{genome}/%s" % config['db']['lastn']['xdir'],
        odir2 = "{genome}/%s" % config['db']['lastp']['xdir'],
        extra = "",
        xpre1 = "{genome}/%s/%s" % (config['db']['lastn']['xdir'], config['db']['lastn']['xpre']),
        xpre2 = "{genome}/%s/%s" % (config['db']['lastp']['xdir'], config['db']['lastp']['xpre']),
        N = "{genome}.%s" % config['last_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['last_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['last_index']['id']),
        j = lambda w: get_resource(w, config, 'last_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'last_index')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        lastdb {params.extra} {params.xpre1} {input.fna}
        lastdb -p {params.extra} {params.xpre2} {input.faa}
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
        j = lambda w: get_resource(w, config, 'bwa_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bwa_index')['ppn']
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
        parallel = lambda w: get_resource(w, config, 'bismark_index')['ppn'] / 2,
        N = "{genome}.%s" % config['bismark_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['bismark_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['bismark_index']['id']),
        j = lambda w: get_resource(w, config, 'bismark_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bismark_index')['ppn']
    conda: "../envs/bismark.yml"
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
        j = lambda w: get_resource(w, config, 'star_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'star_index')['ppn']
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
        j = lambda w: get_resource(w, config, 'gatk_index'),
        mem = lambda w: get_resource(w, config, 'gatk_index')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_index')['ppn']
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
        j = lambda w: get_resource(w, config, 'hisat2_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'hisat2_index')['ppn']
    conda: "../envs/hisat2.yml"
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        hisat2_extract_exons.py {input.gtf} > {params.odir}/db.exon
        hisat2_extract_splice_sites.py {input.gtf} > {params.odir}/db.ss
        hisat2-build -p {threads} --ss {params.odir}/db.ss \
                --exon {params.odir}/db.exon {input.fna} {params.odir}/db
        """

rule salmon_index:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['fna'])
    output:
        "{genome}/%s/%s" % (config['db']['salmon']['xdir'], config['db']['salmon']['xout'])
    params:
        odir = "{genome}/%s" % config['db']['salmon']['xdir'],
        pre = "{genome}/%s/%s" % (config['db']['salmon']['xdir'], config['db']['salmon']['xpre']),
        N = "{genome}.%s" % config['salmon_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['salmon_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['salmon_index']['id']),
        j = lambda w: get_resource(w, config, 'salmon_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'salmon_index')['ppn']
    conda: "../envs/work.yml"
    shell:
        #cut -f1,2 10.tsv | sed '1d' | sort -k1,1 -k2,2 | uniq | awk 'BEGIN{FS="\t";OFS=","}{print $2 $1 $1}' > tx2gene.csv
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        salmon index -p {threads} -t {input.fna} --gencode -i {params.pre}
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
        j = lambda w: get_resource(w, config, 'snpeff_index'),
        mem = lambda w: get_resource(w, config, 'snpeff_index')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'snpeff_index')['ppn']
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

