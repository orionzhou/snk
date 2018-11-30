rule fasta:
    input:
        "{genome}/download/raw.fna"
    output:
        fna = "{genome}/10_genome.fna",
        fai = "{genome}/10_genome.fna.fai",
        chrom_bed = "{genome}/15_intervals/01.chrom.bed",
        chrom_size = "{genome}/15_intervals/01.chrom.sizes",
        gap = "{genome}/15_intervals/11.gap.bed",
    params:
        wdir = lambda w: "%s" % w.genome,
        odir = "08_seq_map",
        prefix = lambda w: config[w.genome]['prefix'],
        N = lambda w: "%s.%s" % (config['fasta']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fasta']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fasta']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'fasta')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'fasta')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'fasta')['mem']
    threads: config['fasta']['ppn']
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
                --merge_short --gap {output.gap} --prefix_chr {params.prefix}

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
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/%s" % (config['blat']['xdir'], config['blat']['x.2bit']),
        "{genome}/21_dbs/%s/%s" % (config['blat']['xdir'], config['blat']['x.ooc']),
    params:
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['blat']['xdir']),
        N = lambda w: "%s.%s" % (config['blat_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['blat_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['blat_index']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'blat_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'blat_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'blat_index')['mem']
    threads: config['blat_index']['ppn']
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        faToTwoBit {input} {output[0]}
        blat {output[0]} tmp.fas tmp.out -makeOoc={output[1]}
        """

rule blast_index:
    input:
        fna = "{genome}/%s/15.fna" % config['adir'],
        faa = "{genome}/%s/15.faa" % config['adir'],
    output:
        "{genome}/21_dbs/%s/%s" % (config['blastn']['xdir'], config['blastn']['xout']),
        "{genome}/21_dbs/%s/%s" % (config['blastp']['xdir'], config['blastp']['xout']),
    params:
        odir1 = lambda w: "%s/21_dbs/%s" % (w.genome, config['blastn']['xdir']),
        odir2 = lambda w: "%s/21_dbs/%s" % (w.genome, config['blastp']['xdir']),
    shell:
        """
        makeblastdb -dbtype nucl -in {input.fna} -title db -out {params.odir1}/db
        makeblastdb -dbtype prot -in {input.faa} -title db -out {params.odir2}/db
        """

rule bwa_index:
    input:
        "{genome}/10_genome.fna"
    output:
        "{genome}/21_dbs/%s/%s" % (config['bwa']['xdir'], config['bwa']['xout'])
    params:
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['bwa']['xdir']),
        N = lambda w: "%s.%s" % (config['bwa_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bwa_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bwa_index']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bwa_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bwa_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bwa_index')['mem']
    threads: config['bwa_index']['ppn']
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        bwa index -a bwtsw -p {params.odir}/db {input}
        """

rule bismark_index:
    input:
        fna = "{genome}/10_genome.fna",
    output:
        "{genome}/21_dbs/%s/%s" % (config['bismark']['xdir'], config['bismark']['xout'])
    params:
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['bismark']['xdir']),
        parallel = lambda w, resources: int(resources.ppn / 2),
        N = lambda w: "%s.%s" % (config['bismark_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bismark_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bismark_index']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bismark_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bismark_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bismark_index')['mem']
    threads: config['bismark_index']['ppn']
    shell:
        """
        rm -rf {params.odir}
        mkdir -p {params.odir}
        cd {params.odir}
        ln -sf ../../10_genome.fna db.fa
        bismark_genome_preparation --bowtie2 --parallel {params.parallel} .
        """

rule star_index:
    input:
        fna = "{genome}/10_genome.fna",
        gtf = "{genome}/%s/10.gtf" % config['adir']
    output:
        "{genome}/21_dbs/%s/%s" % (config['star']['xdir'], config['star']['xout'])
    params:
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['star']['xdir']),
        N = lambda w: "%s.%s" % (config['star_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['star_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['star_index']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'star_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'star_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'star_index')['mem']
    threads: config['star_index']['ppn']
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
        "{genome}/21_dbs/%s/%s" % (config['gatk']['xdir'], config['gatk']['xref']),
        "{genome}/21_dbs/%s/%s" % (config['gatk']['xdir'], config['gatk']['xref.dict']),
    params:
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['gatk']['xdir']),
        N = lambda w: "%s.%s" % (config['gatk_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['gatk_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['gatk_index']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'gatk_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'gatk_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'gatk_index')['mem']
    threads: config['gatk_index']['ppn']
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
        gtf = "{genome}/%s/10.gtf" % config['adir']
    output:
        "{genome}/21_dbs/%s/%s" % (config['hisat2']['xdir'], config['hisat2']['xout'])
    params:
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['hisat2']['xdir']),
        N = lambda w: "%s.%s" % (config['hisat2_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['hisat2_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['hisat2_index']['id'], w.genome),
        q = config['hisat2_index']['q'],
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'hisat2_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'hisat2_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'hisat2_index')['mem']
    threads: config['hisat2_index']['ppn']
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
        fna = "{genome}/10_genome.fna",
        gff = "{genome}/%s/15.gff" % config['adir']
    output:
        "{genome}/21_dbs/%s/%s" % (config['snpeff']['xdir'], config['snpeff']['xcfg']),
        "{genome}/21_dbs/%s/{genome}/%s" % (config['snpeff']['xdir'], config['snpeff']['xout'])
    params:
        fna = lambda w, input: op.abspath(input.fna),
        gff = lambda w, input: op.abspath(input.gff),
        odir = lambda w: "%s/21_dbs/%s" % (w.genome, config['snpeff']['xdir']),
        odir1 = lambda w: "%s/21_dbs/%s/%s" % (w.genome, config['snpeff']['xdir'], w.genome),
        N = lambda w: "%s.%s" % (config['snpeff_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['snpeff_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['snpeff_index']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'snpeff_index')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'snpeff_index')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'snpeff_index')['mem']
    threads: config['snpeff_index']['ppn']
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

