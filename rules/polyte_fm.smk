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
        seqret.py --padding {wildcards.genotype} {input[0]} {output[0]}
        seqret.py --padding {wildcards.genotype} {input[1]} {output[1]}
        """

rule fm2_get_seq_se:
    input:
        "%s/{genotype}_{type}_{opt}.1.fna" % config['fm']['odir2'],
        "%s/{genotype}_{type}_{opt}.2.fna" % config['fm']['odir2'],
    output:
        "%s/{genotype}_{type}_{opt}.se.fna" % config['fm']['odir2'],
    shell:
        """
        fq.merge.py {input[0]} {input[1]} {output}
        """

rule fm2_get_seq_merged:
    input:
        "%s/{genotype}_{type}_{opt}.1.fna" % config['fm']['odir2'],
        "%s/{genotype}_{type}_{opt}.2.fna" % config['fm']['odir2'],
    output:
        "%s/{genotype}_{type}_{opt}.merged.fna" % config['fm']['odir2'],
    shell:
        """
        fq.merge.py --join {input[0]} {input[1]} {output}
        """

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
        sid = lambda wildcards: wildcards.sid,
        pre = lambda wildcards: "%s/%s" % (config['fm']['odir3'], wildcards.sid),
        type = lambda wildcards: config['t'][wildcards.sid]['type'],
        genotype = lambda wildcards: config['t'][wildcards.sid]['genotype'],
        tgt = lambda wildcards: config['t'][wildcards.sid]['tgt'],
        tgt_db = lambda wildcards: "$genome/%s/21_dbs/bwa/db" % config['t'][wildcards.sid]['tgt'],
        opt = lambda wildcards: config['t'][wildcards.sid]['opt'],
        mode = lambda wildcards: config['t'][wildcards.sid]['mode'],
        N = lambda w: "fm3.%s" % (w.sid),
        ppn = config['fm']['bwa']['ppn'],
        walltime = config['fm']['bwa']['walltime'],
        mem = config['fm']['bwa']['mem'],
    threads: config["fm"]['bwa']["ppn"]
    run:
        makedirs(config['fm']['odir3']) 
        if params.mode == 'pe':
            shell("""
            bwa mem -t {threads} -a -T 30 -Y {params.tgt_db} \
                    {input.f1} {input.f2} > {params.pre}.sam
            sam2tsv.py --paired {params.pre}.sam {params.pre}.tsv
            """)
        elif params.mode == 'se':
            shell("""
            bwa mem -t {threads} -a -T 30 -Y {params.tgt_db} \
                    {input.f_se} > {params.pre}.sam
            sam2tsv.py {params.pre}.sam {params.pre}.tsv
            """)
        elif parmas.mode == 'merged':
            shell("""
            bwa mem -t {threads} -a -T 30 -Y {params.tgt_db} \
                    {input.f_merged} > {params.pre}.sam
            sam2tsv.py {params.pre}.sam {params.pre}.tsv
            """)
        shell("""
        atv.filter.py --ident 0.9 --cov 0.9 --best \
                {params.pre}.tsv {params.pre}.filtered.tsv
        """)

rule fm4_coord:
    input:
        "%s/{sid}.sam" % config['fm']['odir3'],
    output:
        "%s/{sid}.tsv" % config['fm']['odir4'],
        "%s/{sid}.filtered.tsv" % config['fm']['odir4']
    params:
        odir = config['fm']['odir4'],
        sid = lambda wildcards: wildcards.sid,
        pre = lambda wildcards: "%s/%s" % (config['fm']['odir4'], wildcards.sid),
        tgt = lambda wildcards: config['t'][wildcards.sid]['tgt'],
        tgt_chain = lambda wildcards: "$genome/%s/08_seq_map/mapb.chain" % config['t'][wildcards.sid]['tgt'],
        N = lambda w: "fm4.%s" % (w.sid),
    shell:
        """
        source activate py27
        mkdir -p {params.odir}
        CrossMap.py bam {params.tgt_chain} {input} - |\
                samtools view -h -O SAM -o {params.pre}.sam
        sam2tsv.py {params.pre}.sam {params.pre}.tsv
        atv.filter.py --ident 0.9 --cov 0.9 --best \
                {params.pre}.tsv {params.pre}.filtered.tsv
        """
 
