rule fm1_get_bed:
    input:
        "%s/{genotype}_{type}.gff" % config['fm']['idir']
    output:
        protected("%s/{genotype}_{type}_{opt}.1.bed" % config['fm']['odir1']),
        protected("%s/{genotype}_{type}_{opt}.2.bed" % config['fm']['odir1']),
    params:
        odir = config['fm']['odir1'],
    shell:
        """
        mkdir -p {params.odir}
        get_flank_bed.R --opt {wildcards.opt} {input} \
                {output[0]} {output[1]}
        """

rule fm2_get_seq:
    input:
        "%s/{genotype}_{type}_{opt}.1.bed" % config['fm']['odir1'],
        "%s/{genotype}_{type}_{opt}.2.bed" % config['fm']['odir1'],
    output:
        "%s/{genotype}_{type}_{opt}.1.fna" % config['fm']['odir2'],
        "%s/{genotype}_{type}_{opt}.2.fna" % config['fm']['odir2'],
    params:
        odir = config['fm']['odir2'],
    shell:
        """
        mkdir -p {params.odir}
        fasta.py extract --padding {wildcards.genotype} {input[0]} >{output[0]}
        fasta.py extract --padding {wildcards.genotype} {input[1]} >{output[1]}
        """

rule fm2_get_seq_se:
    input:
        "%s/{genotype}_{type}_{opt}.1.fna" % config['fm']['odir2'],
        "%s/{genotype}_{type}_{opt}.2.fna" % config['fm']['odir2'],
    output:
        "%s/{genotype}_{type}_{opt}.se.fna" % config['fm']['odir2'],
    shell:
        "fasta.py merge_pe {input[0]} {input[1]} {output}"

rule fm2_get_seq_merged:
    input:
        "%s/{genotype}_{type}_{opt}.1.fna" % config['fm']['odir2'],
        "%s/{genotype}_{type}_{opt}.2.fna" % config['fm']['odir2'],
    output:
        "%s/{genotype}_{type}_{opt}.merged.fna" % config['fm']['odir2'],
    shell:
        "fasta.py merge_pe --join {input[0]} {input[1]} {output}"

def bwa_inputs(wildcards):
    sid = wildcards.sid
    sdic = config['t'][sid]
    genotype, type, opt, mode = sdic['genotype'], sdic['type'], sdic['opt'], sdic['mode']
    pre = "%s/%s_%s_%s" % (config['fm']['odir2'], genotype, type, opt)
    inputs = dict()
    if mode == 'pe':
        inputs['f1'] = "%s.1.fna" % pre
        inputs['f2'] = "%s.2.fna" % pre
    elif mode == 'se':
        inputs['f_se'] = "%s.se.fna" % pre
    elif mode == 'merged':
        inputs['f_se'] = "%s.merged.fna" % pre
    else:
        logging.error("unsupported mode: %s" % mode)
        sys.exit(1)
    return inputs

rule fm3_bwa:
    input:
        unpack(bwa_inputs)
    output:
        "%s/{sid}.sam" % config['fm']['odir3'],
        "%s/{sid}.tsv" % config['fm']['odir3'],
        "%s/{sid}.filtered.tsv" % config['fm']['odir3']
    params:
        odir = config['fm']['odir3'],
        sid = lambda w: w.sid,
        pre = lambda w: "%s/%s" % (config['fm']['odir3'], w.sid),
        type = lambda w: config['t'][w.sid]['type'],
        genotype = lambda w: config['t'][w.sid]['genotype'],
        tgt = lambda w: config['t'][w.sid]['tgt'],
        tgt_db = lambda w: "$genome/%s/21_dbs/bwa/db" % config['t'][w.sid]['tgt'],
        opt = lambda w: config['t'][w.sid]['opt'],
        mode = lambda w: config['t'][w.sid]['mode'],
        N = lambda w: "%s.%s" % (config['fm']['bwa']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fm']['bwa']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fm']['bwa']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'fm', 'bwa')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'fm', 'bwa')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'fm', 'bwa')['mem']
    threads: config["fm"]['bwa']["ppn"]
    run:
        makedirs(config['fm']['odir3'])
        if params.mode == 'pe':
            shell("""
            bwa mem -t {threads} -a -T 30 -Y {params.tgt_db} \
                    {input.f1} {input.f2} > {params.pre}.sam
            sam.py 2tsv --paired {params.pre}.sam >{params.pre}.tsv
            """)
        elif params.mode == 'se':
            shell("""
            bwa mem -t {threads} -a -T 30 -Y {params.tgt_db} \
                    {input.f_se} > {params.pre}.sam
            sam.py 2tsv {params.pre}.sam >{params.pre}.tsv
            """)
        elif parmas.mode == 'merged':
            shell("""
            bwa mem -t {threads} -a -T 30 -Y {params.tgt_db} \
                    {input.f_merged} > {params.pre}.sam
            sam.py 2tsv {params.pre}.sam >{params.pre}.tsv
            """)
        shell("""
        atv.py filter --ident 0.9 --cov 0.9 --best \
                {params.pre}.tsv >{params.pre}.filtered.tsv
        """)

rule fm4_coord:
    input:
        "%s/{sid}.sam" % config['fm']['odir3'],
    output:
        "%s/{sid}.tsv" % config['fm']['odir4'],
        "%s/{sid}.filtered.tsv" % config['fm']['odir4']
    params:
        odir = config['fm']['odir4'],
        sid = lambda w: w.sid,
        pre = lambda w: "%s/%s" % (config['fm']['odir4'], w.sid),
        tgt = lambda w: config['t'][w.sid]['tgt'],
        tgt_chain = lambda w: "$genome/%s/08_seq_map/mapb.chain" % config['t'][w.sid]['tgt'],
        N = lambda w: "%s.%s" % (config['fm']['coord']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['fm']['coord']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['fm']['coord']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'fm', 'coord')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'fm', 'coord')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'fm', 'coord')['mem']
    shell:
#        CrossMap.py bam {params.tgt_chain} {input} - |\
#                samtools view -h -O SAM -o {params.pre}.sam
#        source activate py27
#        source deactivate
        """
        mkdir -p {params.odir}

        sam.py 2tsv --paired {params.pre}.sam >{params.pre}.tsv
        atv.py filter --ident 0.9 --cov 0.9 --best \
                {params.pre}.tsv >{params.pre}.filtered.tsv
        """

