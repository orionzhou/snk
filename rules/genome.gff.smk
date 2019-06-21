def fix_option(genome):
    if genome in 'Mo17 W22 PH207 PHB47 HZS'.split():
        return genome.lower()
    else:
        return 'ensembl'

def anno1_input(w):
    genome = w.genome
    if config['x'][genome]['hybrid']:
        genome1, genome2 = genome.split('x')
        return dict(
            gff1 = "%s/%s/%s" % (genome1, config['db']['annotation']['xdir'], config['db']['annotation']['gff']),
            gff2 = "%s/%s/%s" % (genome2, config['db']['annotation']['xdir'], config['db']['annotation']['gff'])
        )
    else:
        return dict(
            chain = "%s/%s/%s" % (genome, config['db']['fasta']['xdir'], config['db']['fasta']['fchain']),
            gff = "{genome}/download/raw.gff"
        )

rule anno1_clean:
    input: unpack(anno1_input)
    output:
        "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['gff']),
    params:
        gdir = "{genome}",
        wdir = "{genome}/%s" % config['db']['annotation']['xdir'],
        fixopt = lambda w: fix_option(w.genome),
        N = "%s.{genome}" % config['anno1_clean']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno1_clean']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno1_clean']['id']),
        j = lambda w: get_resource(w, config, 'anno1_clean'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'anno1_clean')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/anno1.clean.py"

rule anno2_index:
    input:
        ref = "{genome}/%s/%s" % (config['db']['fasta']['xdir'], config['db']['fasta']['ref']),
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
        N = "%s.{genome}" % config['anno2_index']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno2_index']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno2_index']['id']),
        j = lambda w: get_resource(w, config, 'anno2_index'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'anno2_index')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        gff.py index {input.gff} {output.db}
        gff.py 2tsv {input.gff} > {output.tsv}
        gff.py 2gtf {input.gff} > {output.gtf}
        gff.py 2bed12 {input.gff} > {output.bed}
        gff.py note --attribute note1,note2 {input.gff} > {output.des}
        gff.py 2fas {input.gff} {input.ref} >{output.fna}
        fasta.py translate {output.fna} > {output.faa}

        gff.py picklong {input.gff} > {output.pgff}
        gff.py index {output.pgff} {output.pdb}
        gff.py 2tsv {output.pgff} > {output.ptsv}
        gff.py 2gtf {output.pgff} > {output.pgtf}
        gff.py 2bed12 {output.pgff} > {output.pbed}
        gff.py note --attribute note1,note2 {output.pgff} > {output.pdes}
        gff.py 2fas {output.pgff} {input.ref} > {output.pfna}
        fasta.py translate {output.pfna} > {output.pfaa}
        """

rule anno3_blast:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
        blastn_db = "{genome}/%s/%s" % (config['db']['blastn']['xdir'], config['db']['blastn']['xout']),
        blastp_db = "{genome}/%s/%s" % (config['db']['blastp']['xdir'], config['db']['blastp']['xout']),
    output:
        cds = "{genome}/%s/01.blastn.tsv" % config['db']['tandup']['xdir'],
        pro = "{genome}/%s/01.blastp.tsv" % config['db']['tandup']['xdir'],
    params:
        gdir = "{genome}",
        blastn_db = "{genome}/%s/%s" % (config['db']['blastn']['xdir'], config['db']['blastn']['xpre']),
        blastp_db = "{genome}/%s/%s" % (config['db']['blastp']['xdir'], config['db']['blastp']['xpre']),
        N = "%s.{genome}" % config['anno3_blast']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno3_blast']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno3_blast']['id']),
        j = lambda w: get_resource(w, config, 'anno3_blast'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'anno3_blast')['ppn']
    conda: "../envs/blast.yml"
    shell:
        """
        blastn -db {params.blastn_db} -query {input.fna} -num_threads {threads} \
            -outfmt 6 -out {output.cds}
        blastp -db {params.blastp_db} -query {input.faa} -num_threads {threads} \
            -outfmt 6 -out {output.pro}
        """

rule anno4_tandup:
    input:
        fna = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfna']),
        faa = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lfaa']),
        gbed = "{genome}/%s/%s" % (config['db']['annotation']['xdir'], config['db']['annotation']['lbed']),
        cds = "{genome}/%s/01.blastn.tsv" % config['db']['tandup']['xdir'],
        pro = "{genome}/%s/01.blastp.tsv" % config['db']['tandup']['xdir'],
    output:
        cds = "{genome}/%s/%s" % (config['db']['tandup']['xdir'], config['db']['tandup']['xout'].replace('pro','cds')),
        pro = "{genome}/%s/%s" % (config['db']['tandup']['xdir'], config['db']['tandup']['xout']),
    params:
        N = "%s.{genome}" % config['anno4_tandup']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['anno4_tandup']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['anno4_tandup']['id']),
        j = lambda w: get_resource(w, config, 'anno4_tandup'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'anno4_tandup')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        python -m jcvi.compara.catalog tandem --strip_gene_name=None {input.cds} \
            {input.fna} {input.gbed} > {output.cds}
        python -m jcvi.compara.catalog tandem --strip_gene_name=None {input.pro} \
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
        "{genome}/%s/%s" % (config['db']['rds']['xdir'], config['db']['rds']['xout'])
    params:
        wdir = "{genome}",
        N = "%s.{genome}" % config['prepR']['id'],
        e = "{genome}/%s/%s.e" % (config['dirj'], config['prepR']['id']),
        o = "{genome}/%s/%s.o" % (config['dirj'], config['prepR']['id']),
        j = lambda w: get_resource(w, config, 'prepR'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'prepR')['ppn']
    conda: "../envs/work.yml"
    shell: "genome.prep.R {wildcards.genome}"



