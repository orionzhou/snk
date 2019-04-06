def fix_option(genome):
    if genome in 'Mo17 W22 PH207 PHB47'.split():
        return genome.lower()
    else:
        return 'ensembl'

rule anno1_clean:
    input:
        chain = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['fchain']),
        gff = "{genome}/download/raw.gff"
    output:
        "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gff']),
    params:
        gdir = "{genome}",
        wdir = "{genome}/%s" % config['db']['annotation']['xdir'],
        prefix = lambda w: config['x'][w.genome]['prefix'],
        fixopt = lambda w: fix_option(w.genome),
        N = "%s.{genome}" % config['anno1_clean']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno1_clean']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno1_clean']['id']),
    conda: "../envs/work.yml"
    shell:
        """
        gff.py fix --opt {params.fixopt} {input.gff} > {params.wdir}/01.fixed.gff
        liftOver -gff {params.wdir}/01.fixed.gff {input.chain} \
                {params.wdir}/02.lifted.gff {params.wdir}/unmapped
        ln -sf 02.lifted.gff {output}
        """

rule anno2_index:
    input:
        gff = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gff']),
    output:
        db = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gff_db']),
        gtf = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gtf']),
        tsv = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['tsv']),
        bed = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['bed']),
        des = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['des']),
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['fna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['faa']),
        pgff = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lgff']),
        pdb = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lgff_db']),
        pgtf = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lgtf']),
        ptsv = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['ltsv']),
        pbed = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lbed']),
        pdes = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['ldes']),
        pfna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        pfaa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
    params:
        gdir = "{genome}",
        wdir = "{genome}/%s" % config['adir'],
        db = "{genome}/%s/10.gff.db" % config['adir'],
        ref = "{genome}/10_genome.fna",
        N = "%s.{genome}" % config['anno2_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno2_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno2_index']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'anno2_index')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'anno2_index')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'anno2_index')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'anno2_index')['mem']
    threads: config['anno2_index']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        gff.py index {input.gff} {output.db}
        gff.py 2tsv {input.gff} > {output.tsv}
        gff.py 2gtf {input.gff} > {output.gtf}
        gff.py 2bed12 {input.gff} > {output.bed}
        gff.py note --attribute note1,note2 {input.gff} > {output.des}
        gff.py 2fas {input.gff} {params.ref} >{output.fna}
        fasta.py translate {output.fna} > {output.faa}

        gff.py picklong {input.gff} > {output.pgff}
        gff.py index {output.pgff} {output.pdb}
        gff.py 2tsv {output.pgff} > {output.ptsv}
        gff.py 2gtf {output.pgff} > {output.pgtf}
        gff.py 2bed12 {output.pgff} > {output.pbed}
        gff.py note --attribute note1,note2 {output.pgff} > {output.pdes}
        gff.py 2fas {output.pgff} {params.ref} > {output.pfna}
        fasta.py translate {output.pfna} > {output.pfaa}
        """

rule anno3_blast:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
        blastn_db = "{genome}/%s/%s" % (config['db']['blastn']['xdir'], config['db']['blastn']['xout']),
        blastp_db = "{genome}/%s/%s" % (config['db']['blastp']['xdir'], config['db']['blastp']['xout']),
    output:
        cds = "{genome}/%s/21.blast.cds.tsv" % config['db']['annotation']['xdir'],
        pro = "{genome}/%s/21.blast.pro.tsv" % config['db']['annotation']['xdir'],
    params:
        gdir = "{genome}",
        wdir = "{genome}/%s" % config['db']['annotation']['xdir'],
        N = "%s.{genome}" % config['anno3_blast']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno3_blast']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno3_blast']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'anno3_blast')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'anno3_blast')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'anno3_blast')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'anno3_blast')['mem']
    threads: config['anno3_blast']['ppn']
    conda: "../envs/blast.yml"
    shell:
        """
        blastn -db {input.blastn_db} -query {input.fna} -num_threads {threads} \
            -outfmt 6 -out {output.cds}
        blastp -db {input.blastp_db} -query {input.faa} -num_threads {threads} \
            -outfmt 6 -out {output.pro}
        """

rule anno4_tandup:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
        gbed = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lbed']),
        cds = "{genome}/%s/21.blast.cds.tsv" % config['db']['annotation']['xdir'],
        pro = "{genome}/%s/21.blast.pro.tsv" % config['db']['annotation']['xdir'],
    output:
        cds = "{genome}/%s/22.tandup.cds.tsv" % config['db']['annotation']['xdir'],
        pro = "{genome}/%s/22.tandup.pro.tsv" % config['db']['annotation']['xdir'],
    params:
        gdir = "{genome}",
        wdir = "{genome}/%s" % config['adir'],
        N = "%s.{genome}" % config['anno4_tandup']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno4_tandup']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno4_tandup']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'anno4_tandup')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'anno4_tandup')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'anno4_tandup')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'anno4_tandup')['mem']
    threads: config['anno4_tandup']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        catalog.py tandem --strip_gene_name None {input.cds} \
            {input.fna} {input.gbed} > {output.cds}
        catalog.py tandem --strip_gene_name None {input.pro} \
            {input.faa} {input.gbed} > {output.pro}
        """

rule prepR:
    input:
        chrom_size = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['chrom_size']),
        chrom_bed = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['chrom_bed']),
        gap_bed = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['gap']),
        gene_tsv = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['ltsv']),
        gene_des = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['ldes']),
    output:
        "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['rds'])
    params:
        wdir = "{genome}",
        N = "%s.{genome}" % config['prepR']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['prepR']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['prepR']['id']),
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'prepR')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'prepR')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'prepR')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'prepR')['mem']
    threads: config['prepR']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        genome.prep.R {wildcards.genome}
        """
