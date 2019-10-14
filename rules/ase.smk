rule ase3_snp:
    input:
        bam = "{yid}/%s/{sid}.bam" % config['rnaseq']['idir'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]["fasta"]['ref'],
    output:
        "{yid}/%s/{sid}.tsv.gz" % config['rnaseq']['od33c']
    params:
        bcf = lambda w: "%s/%s.bcf" % (config['rnaseq']['ase_bcf'], config['y'][w.yid]['t'][w.sid]['Genotype']),
        mq = 10,
        N = "{yid}.%s.{sid}" % config['ase3_snp']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['ase3_snp']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['ase3_snp']['id']),
        j = lambda w: get_resource(w, config, 'ase3_snp'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase3_snp')['ppn']
    conda: "../envs/alfred.yml"
    shell:
        """
        alfred ase -r {input.ref} -s sample -v {params.bcf} \
            -m {params.mq} -p -a {output} {input.bam}
        """

rule ase4_split:
    input:
        bam = "{yid}/%s/{sid}.bam" % config['rnaseq']['idir'],
#        bam = "{yid}/%s/{sid}.bam" % config['rnaseq']['od33c'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]["fasta"]['ref'],
    output:
        bam1 = "{yid}/%s/{sid}.1.bam" % config['rnaseq']['od33d'],
        bam2 = "{yid}/%s/{sid}.2.bam" % config['rnaseq']['od33d'],
        tsv1 = "{yid}/%s/{sid}.1.tsv" % config['rnaseq']['od33d'],
        tsv2 = "{yid}/%s/{sid}.2.tsv" % config['rnaseq']['od33d']
    params:
        bcf = lambda w: "%s/%s.bcf" % (config['rnaseq']['ase_bcf'], config['y'][w.yid]['t'][w.sid]['Genotype']),
        mq = 10,
        gtf = lambda w: db_index(w, 'gtf'),
        extra = featurecounts_extra,
        strand = mmquant_strand,
        N = "{yid}.%s.{sid}" % config['ase4_split']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['ase4_split']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['ase4_split']['id']),
        j = lambda w: get_resource(w, config, 'ase4_split'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase4_split')['ppn']
    conda: "../envs/alfred.yml"
    shell:
#        mmquant -a {params.gtf} -r {output.bam1} \
#            -o {output.tsv1} -l 0.3 -t {threads} -s {params.strand}
#        mmquant -a {params.gtf} -r {output.bam2} \
#            -o {output.tsv2} -l 0.3 -t {threads} -s {params.strand}
        """
        alfred split -r {input.ref} -s sample -v {params.bcf} \
            -m {params.mq} -p {output.bam1} -q {output.bam2} {input.bam}
        featureCounts -T {threads} {params.extra} \
            -a {params.gtf} -o {output.tsv1} {output.bam1}
        featureCounts -T {threads} {params.extra} \
            -a {params.gtf} -o {output.tsv2} {output.bam2}
        """

def ase5_inputs(w):
    yid = w.yid
    cnts, cnts2 = [], []
    for sid in config['y'][yid]['SampleID']:
        gt = config['y'][yid]['t'][sid]['Genotype']
        ase_bcf = "%s/%s.bcf" % (config['rnaseq']['ase_bcf'], gt)
        if op.isfile(ase_bcf):
            cnts.append("%s/%s/%s.%s.tsv" % (yid, config['rnaseq']['od33d'], sid, '1'))
            cnts.append("%s/%s/%s.%s.tsv" % (yid, config['rnaseq']['od33d'], sid, '2'))
            cnts2.append("%s/%s/%s.tsv.gz" % (yid, config['rnaseq']['od33c'], sid))
        else:
            print("  %s|%s [%s] skipped" % (yid, sid, gt))
#    stats = expand("%s/%s/{sid}.tsv" % (yid, config['rnaseq']['od33c']), sid = config['y'][yid]['SampleID'])
#    cnts = expand("%s/%s/{sid}.{suf}.tsv" % (yid, config['rnaseq']['od33d']), sid = config['y'][yid]['SampleID'], suf = ['1','2'])
    return dict(cnts=cnts, cnts2=cnts2)

rule ase5_merge:
    input: unpack(ase5_inputs)
    output:
        cnt = protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_ase'])),
        cnt2 = protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_ase2']))
    params:
        N = "{yid}.%s" % config['ase5_merge']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['ase5_merge']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['ase5_merge']['id']),
        j = lambda w: get_resource(w, config, 'ase5_merge'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase5_merge')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        merge.stats.R --opt ase -o {output.cnt} {input.cnts}
        merge.stats.R --opt ase_snp -o {output.cnt2} {input.cnts2}
        """

