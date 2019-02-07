def fix_option(genome):
    if genome in 'Mo17 W22 PH207 PHB47'.split():
        return genome.lower()
    else:
        return 'ensembl'

rule anno1_clean:
    input:
        chain = "{genome}/08_seq_map/mapf.chain",
        gff = "{genome}/download/raw.gff"
    output:
        "{genome}/%s/10.gff" % config['adir']
    params:
        gdir = lambda w: "%s" % w.genome,
        wdir = lambda w: "%s/%s" % (w.genome, config['adir']),
        prefix = lambda w: config[w.genome]['prefix'],
        fixopt = lambda w: fix_option(w.genome),
    shell:
        """
        gff.py fix --opt {params.fixopt} {input.gff} > {params.wdir}/01.fixed.gff
        liftOver -gff {params.wdir}/01.fixed.gff {input.chain} \
                {params.wdir}/02.lifted.gff {params.wdir}/unmapped
        ln -sf 02.lifted.gff {output}
        """

rule anno2_index:
    input:
        gff = "{genome}/%s/10.gff" % config['adir']
    output:
        db = "{genome}/%s/10.gff.db" % config['adir'],
        tsv = "{genome}/%s/10.tsv" % config['adir'],
        gtf = "{genome}/%s/10.gtf" % config['adir'],
        bed = "{genome}/%s/10.bed" % config['adir'],
        des = "{genome}/%s/10.desc.tsv" % config['adir'],
        fna = "{genome}/%s/10.fna" % config['adir'],
        faa = "{genome}/%s/10.faa" % config['adir'],
        pgff = "{genome}/%s/15.gff" % config['adir'],
        pdb = "{genome}/%s/15.gff.db" % config['adir'],
        ptsv = "{genome}/%s/15.tsv" % config['adir'],
        pgtf = "{genome}/%s/15.gtf" % config['adir'],
        pbed = "{genome}/%s/15.bed" % config['adir'],
        pdes = "{genome}/%s/15.desc.tsv" % config['adir'],
        pfna = "{genome}/%s/15.fna" % config['adir'],
        pfaa = "{genome}/%s/15.faa" % config['adir'],
    params:
        gdir = lambda w: "%s" % w.genome,
        wdir = lambda w: "%s/%s" % (w.genome, config['adir']),
        db = lambda w: "%s/%s/10.gff.db" % (w.genome, config['adir']),
        ref = lambda w: "%s/10_genome.fna" % w.genome,
        N = lambda w: "%s.%s" % (config['anno2_index']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['anno2_index']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['anno2_index']['id'], w.genome),
        q = lambda w, resources: resources.q,
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'anno2_index')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'anno2_index')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'anno2_index')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'anno2_index')['mem']
    threads: config['anno2_index']['ppn']
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

rule anno3_tandup:
    input:
        fna = "{genome}/%s/15.fna" % config['adir'],
        faa = "{genome}/%s/15.faa" % config['adir'],
        gbed = "{genome}/%s/15.bed" % config['adir'],
        blastn_db = "{genome}/21_dbs/%s/%s" % (config['blastn']['xdir'], config['blastn']['xout']),
        blastp_db = "{genome}/21_dbs/%s/%s" % (config['blastp']['xdir'], config['blastp']['xout']),
    output:
        cds = "{genome}/%s/22.tandup.cds.tsv" % config['adir'],
        pro = "{genome}/%s/22.tandup.pro.tsv" % config['adir'],
    params:
        gdir = lambda w: "%s" % w.genome,
        wdir = lambda w: "%s/%s" % (w.genome, config['adir']),
        ref = lambda w: "%s/10_genome.fna" % w.genome,
        blastn_db = lambda w: "%s/21_dbs/%s/%s" % (w.genome, config['blastn']['xdir'], config['blastn']['xpre']),
        blastp_db = lambda w: "%s/21_dbs/%s/%s" % (w.genome, config['blastp']['xdir'], config['blastp']['xpre']),
        N = lambda w: "%s.%s" % (config['anno3_tandup']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['anno3_tandup']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['anno3_tandup']['id'], w.genome),
        q = lambda w, resources: resources.q,
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'anno3_tandup')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'anno3_tandup')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'anno3_tandup')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'anno3_tandup')['mem']
    threads: config['anno3_tandup']['ppn']
    shell:
        """
        blastn -db {params.blastn_db} -query {input.fna} -num_threads {threads} \
            -outfmt 6 -out {params.wdir}/21.blast.cds.tsv
        blastp -db {params.blastp_db} -query {input.faa} -num_threads {threads} \
            -outfmt 6 -out {params.wdir}/21.blast.pro.tsv
        catalog.py tandem --strip_gene_name None {params.wdir}/21.blast.cds.tsv \
            {input.fna} {input.gbed} > {output.cds}
        catalog.py tandem --strip_gene_name None {params.wdir}/21.blast.pro.tsv \
            {input.faa} {input.gbed} > {output.pro}
        """
