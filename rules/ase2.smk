def ase1_extra(w):
    yid, sid = w.yid, w.sid
    extra = "--is_sorted"
    dir_h5 = '/home/springer/zhoux379/projects/reseq/ase/B73xMo17/07_h5'
    extra += " --snp_tab %s/snp_tab.h5" % dir_h5
    extra += " --snp_index %s/snp_index.h5" % dir_h5
    extra += " --haplotype %s/haplotypes.h5" % dir_h5
    if config['y'][yid]['t'][sid]['paired']:
        extra += " --is_paired_end"
    return extra

def ase1_output_str(w):
    yid, sid = w.yid, w.sid
    pre = "%s/%s/%s.remap" % (yid, config['rnaseq']['od33a'], sid)
    fq = "%s.fq.gz" % pre
    fq1, fq2 = ("%s.fq%s.gz" % (pre, suf) for suf in '1 2'.split())
    output_str = ''
    if config['y'][yid]['t'][sid]['paired']:
        output_str = "-fq %s -fq2 %s" % (fq1, fq2)
    else:
        output_str = "-fq %s" % fq
    return output_str

def ase1_cmds(w):
    yid, sid = w.yid, w.sid
    cmds = []
    odir = "%s/%s" % (yid, config['rnaseq']['od33a'])
    if config['y'][yid]['t'][sid]['paired']:
        cmds.append("zcat %s/%s.remap.fq1.gz | wc -l" % (odir, sid))
        cmds.append("zcat %s/%s.remap.fq2.gz | wc -l" % (odir, sid))
    else:
        cmds.append("zcat %s/%s.remap.fq.gz | wc -l" % (odir, sid))
    return "\n".join(cmds)

rule ase1_find_intersect:
    input: "{yid}/%s/{sid}.bam" % config['rnaseq']['idir']
    output:
        keep = "{yid}/%s/{sid}.keep.bam" % config['rnaseq']['od33a'],
        remap = "{yid}/%s/{sid}.to.remap.bam" % config['rnaseq']['od33a'],
    params:
        odir = "{yid}/%s" % config['rnaseq']['od33a'],
        extra = ase1_extra,
        cmd = ase1_cmds,
        qsort = "{yid}/%s/{sid}.to.remap.qsort.bam" % config['rnaseq']['od33a'],
        output_str = ase1_output_str,
        N = "{yid}.%s.{sid}" % config['ase1_find_intersect']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['ase1_find_intersect']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['ase1_find_intersect']['id']),
        j = lambda w: get_resource(w, config, 'ase1_find_intersect'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase1_find_intersect')['ppn']
    conda: "../envs/wasp.yml"
    shell:
#        rm {params.odir}/{wildcards.sid}.*.gz
#        samtools sort -n -o {params.qsort} {output.remap}
#        bedtools bamtofastq -i {params.qsort} {params.output_str}
        """
        python $src/git/WASP/mapping/find_intersecting_snps.py \
            {params.extra} \
            --output_dir {params.odir} \
            {input}
        {params.cmd}
        """

def ase2_input_str(w):
    yid, sid = w.yid, w.sid
    pre = "%s/%s/%s.remap" % (yid, config['rnaseq']['od33a'], sid)
    fq = "%s.fq.gz" % pre
    fq1, fq2 = ("%s.fq%s.gz" % (pre, suf) for suf in '1 2'.split())
    input_str = ''
    if config['y'][yid]['t'][sid]['paired']:
        input_str = "-1 %s -2 %s" % (fq1, fq2)
    else:
        input_str = "-U %s" % fq
    return input_str

rule ase2_remap:
    input: "{yid}/%s/{sid}.to.remap.bam" % config['rnaseq']['od33a']
    output:
        temp("{yid}/%s/{sid}.bam" % config['rnaseq']['od33b']),
        "{yid}/%s/{sid}.txt" % config['rnaseq']['od33b'],
        "{yid}/%s/{sid}.sorted.bam" % config['rnaseq']['od33b']
    params:
        input_str = ase2_input_str,
        index = lambda w: db_index(w, 'hisat2'),
        extra = hisat2_extra,
        tmpdir = config['tmpdir'],
        N = "{yid}.%s.{sid}" % config['ase2_remap']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['ase2_remap']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['ase2_remap']['id']),
        j = lambda w: get_resource(w, config, 'ase2_remap'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase2_remap')['ppn']
    conda: "../envs/hisat2.yml"
    shell:
        """
        hisat2 {params.extra} --threads {threads} \
            -x {params.index} {params.input_str} \
            --summary-file {output[1]} \
            | samtools view -Sbh -o {output[0]} -
        sambamba sort -t {threads} --tmpdir={params.tmpdir} \
            -o {output[2]} {output[0]}
        """

rule ase3_filter:
    input:
        keep = "{yid}/%s/{sid}.keep.bam" % config['rnaseq']['od33a'],
        before = "{yid}/%s/{sid}.to.remap.bam" % config['rnaseq']['od33a'],
        after = "{yid}/%s/{sid}.sorted.bam" % config['rnaseq']['od33b'],
    output:
        bam = "{yid}/%s/{sid}.bam" % config['rnaseq']['od33c'],
        tsv = "{yid}/%s/{sid}.tsv" % config['rnaseq']['od33c']
    params:
        keep = "{yid}/%s/{sid}.keep.bam" % config['rnaseq']['od33c'],
        unsort = "{yid}/%s/{sid}.unsorted.bam" % config['rnaseq']['od33c'],
        tmpdir = config['tmpdir'],
        N = "{yid}.%s.{sid}" % config['ase3_filter']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['ase3_filter']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['ase3_filter']['id']),
        j = lambda w: get_resource(w, config, 'ase3_filter'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ase3_filter')['ppn']
    conda: "../envs/wasp.yml"
    shell:
        """
        python $src/git/WASP/mapping/filter_remapped_reads.py \
            {input.before} {input.after} {params.keep}
        samtools merge -f \
            {params.unsort} {params.keep} {input.keep}
        sambamba sort -t {threads} --tmpdir={params.tmpdir} \
            -o {output.bam} {params.unsort}
        bam.py stat {output.bam} > {output.tsv}
        """


