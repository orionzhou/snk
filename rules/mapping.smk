def mapping_inputs(w):
    yid, sid = w.yid, w.sid
    part = w.part
    pre = "%s_%s" % (sid, part)
    idir = config['mapping']['idir']

    inputs = dict()
    if config['y'][yid]['t'][sid]['paired']:
        inputs['fq1'] = ancient("%s/%s/%s_1.fq.gz" % (yid, idir, pre))
        inputs['fq2'] = ancient("%s/%s/%s_2.fq.gz" % (yid, idir, pre))
    else:
        inputs['fq'] = ancient("%s/%s/%s.fq.gz" % (yid, idir, pre))
    return inputs

def mapping_input_str(w, input, mapper):
    yid, sid = w.yid, w.sid
    assert mapper in config['valid']['mapper'], "invalid mapper: %s" % mapper
    input_str = ''
    if config['y'][yid]['t'][sid]['paired']:
        if mapper in ['hisat2','bismark']:
            input_str = "-1 %s -2 %s" % (input['fq1'], input['fq2'])
        elif mapper in ['star','bwa']:
            input_str = "%s %s" % (input['fq1'], input['fq2'])
    else:
        if mapper in ['hisat2']:
            input_str = "-U %s" % input['fq']
        else:
            input_str = "%s" % input['fq']
    return input_str

def db_index(w, db):
    yid, sid = w.yid, w.sid
    ref = config['y'][yid]['ref']
#    if 'ase' in config['y'][yid] and config['y'][yid]['ase']:
#        ref = config['y'][yid]['t'][sid]['Genotype']
#        ref = normalize_genotype(ref)

    if db == 'star':
        return config['g'][ref]["star"]['xpre']
    elif db == 'hisat2':
        return config['g'][ref]["hisat2"]['xpre']
    elif db == 'bwa':
        return config['g'][ref]["bwa"]['xpre']
    elif db == 'bismark':
        return config['g'][ref]["bismark"]['xpre']
    elif db == 'gtf':
        return config['g'][ref]["annotation"]['gtf']
    else:
        print('unknown db: %d' % db)
        sys.exit(1)

def star_extra(w):
    extras = """
        --outSAMmapqUnique 60
        --outFilterType BySJout
        --outFilterMultimapNmax 20
        --alignSJoverhangMin 8
        --alignSJDBoverhangMin 1
        --outFilterMismatchNmax 999
        --outFilterMismatchNoverReadLmax 1.0
        --alignIntronMin 20
        --alignIntronMax 1000000
        --alignMatesGapMax 1000000
        --outSAMunmapped Within KeepPairs
    """.split()
    extras.append("--outSAMattrRGline ID:%s SM:%s" % (w.sid, w.sid))
    #if 'vcf' in config and wildcards.gt in config['vcf']:
    if 1 == 2:
        extras.append("--varVCFfile %s" % config['vcf'][w.gt])
        extras.append("--waspOutputMode SAMtag")
        extras.append("--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch vG vA vW")
    else:
        extras.append("--outSAMattributes All")
    return " ".join(extras)

rule star:
    input: unpack(mapping_inputs)
    output:
        temp("{yid}/%s/{sid}_{part}/Aligned.out.bam" % config['mapping']['od21s']),
        "{yid}/%s/{sid}_{part}/Log.final.out" % config['mapping']['od21s']
    params:
        index = lambda w: db_index(w, 'star'),
        input_str = lambda w, input: mapping_input_str(w, input, 'star'),
        outprefix = "{yid}/%s/{sid}_{part}/" % config['mapping']['od21s'],
        readcmd = "--readFilesCommand zcat",
        extra = star_extra,
        N = "{yid}.%s.{sid}{part}" % config['star']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['star']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['star']['id']),
        j = lambda w: get_resource(w, config, 'star'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'star')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix} \
        --outStd Log
        """

def hisat2_extra(w):
    extras = "--new-summary --no-softclip".split()
    extras = "--new-summary".split()
    extras.append("--rg-id %s --rg SM:%s" % (w.sid, w.sid))
    return " ".join(extras)

rule hisat2:
    input: unpack(mapping_inputs)
    output:
        temp("{yid}/%s/{sid}_{part}.bam" % config['mapping']['od21h']),
        "{yid}/%s/{sid}_{part}.txt" % config['mapping']['od21h']
    params:
        index = lambda w: db_index(w, 'hisat2'),
        input_str = lambda w, input: mapping_input_str(w, input, 'hisat2'),
        extra = hisat2_extra,
        N = "{yid}.%s.{sid}{part}" % config['hisat2']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['hisat2']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['hisat2']['id']),
        j = lambda w: get_resource(w, config, 'hisat2'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'hisat2')['ppn']
    conda: "../envs/hisat2.yml"
    shell:
        """
        hisat2 {params.extra} --threads {threads} \
            -x {params.index} {params.input_str} \
            --summary-file {output[1]} \
            | samtools view -Sbh -o {output[0]} -
        """

def bwa_extra(w):
    extras = []
    yid, sid = w.yid, w.sid
    sm = sid
    if 'Genotype' in config['y'][yid]['t'][sid]:
        sm = config['y'][yid]['t'][w.sid]['Genotype']
    pl = 'ILLUMINA'
    if config['y'][yid]['readtype'] == 'solid':
        pl = 'SOLID'
    extras.append("-R '@RG\\tID:%s\\tSM:%s\\tPL:%s'" % (sid, sm, pl))
    return " ".join(extras)

rule bwa:
    input: unpack(mapping_inputs)
    output: temp("{yid}/%s/{sid}_{part}.bam" % config['mapping']['od21b'])
    params:
        index = lambda w: db_index(w, 'bwa'),
        input_str = lambda w, input: mapping_input_str(w, input, 'bwa'),
        extra = bwa_extra,
        N = "{yid}.%s.{sid}{part}" % config['bwa']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['bwa']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['bwa']['id']),
        j = lambda w: get_resource(w, config, 'bwa'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bwa')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        bwa mem -t {threads} {params.index} {params.extra} \
            {params.input_str} |\
            samtools view -bS - > {output}
        """

def bismark_extra(w):
    yid, sid = w.yid, w.sid
    sm = sid
    extras = ["-n 1"]
    if 'Genotype' in config['y'][yid]['t'][sid]:
        sm = config['y'][yid]['t'][w.sid]['Genotype']
    extras.append("--rg_tag --rg_id %s --rg_sample %s" % (sid, sm))
    return " ".join(extras)

rule bismark:
    input: unpack(mapping_inputs)
    output:
        bam = "{yid}/%s/{sid}_{part}.bam" % config['mapping']['od21m'],
        report = "{yid}/%s/{sid}_{part}.txt" % config['mapping']['od21m'],
    params:
        index = lambda w: db_index(w, 'bismark'),
        input_str = lambda w, input: mapping_input_str(w, input, 'bismark'),
        odir = "{yid}/%s" % config['mapping']['od21m'],
        opre = '{sid}_{part}',
        tmp_dir = "{yid}/%s/{sid}_{part}" % config['mapping']['od21m'],
        p = lambda w: int(get_resource(w, config, 'bismark')['ppn'] / 2),
        extra = bismark_extra,
        N = "{yid}.%s.{sid}{part}" % config['bismark']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['bismark']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['bismark']['id']),
        j = lambda w: get_resource(w, config, 'bismark'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bismark')['ppn']
    conda: "../envs/bismark.yml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        bismark --parallel {params.p} \
            {params.index} {params.input_str} \
            --temp_dir {params.tmp_dir} \
            {params.extra} --output_dir {params.odir}
        mv {params.odir}/{params.opre}*bismark*.bam {output.bam}
        mv {params.odir}/{params.opre}*bismark*_report.txt {output.report}
        """

def sambamba_sort_input(w):
    yid, sid = w.yid, w.sid
    parts = config['y'][yid]['t'][sid]['parts']
    mapper = config['y'][yid]['mapper']
    assert mapper in config['valid']['mapper'], "invalid mapper: %s" % mapper
    ptn = ''
    if mapper == 'star':
        ptn = "%s/%s/%s_{part}/Aligned.out.bam" % (yid, config['mapping']['od21s'], sid)
    elif mapper == 'hisat2':
        ptn = "%s/%s/%s_{part}.bam" % (yid, config['mapping']['od21h'], sid)
    elif mapper == 'bismark':
        ptn = "%s/%s/%s_{part}.bam" % (yid, config['mapping']['od21m'], sid)
    elif mapper == 'bwa':
        ptn = "%s/%s/%s_{part}.bam" % (yid, config['mapping']['od21b'], sid)
    fis = expand(ptn, part = parts)
    return fis

rule sambamba_sort:
    input: sambamba_sort_input
    output:
        "{yid}/%s/{sid}.bam" % config['mapping']['odir'],
        "{yid}/%s/{sid}.bam.bai" % config['mapping']['odir']
    params:
        mapper = lambda w: config['y'][w.yid]['mapper'],
        odir = "{yid}/%s" % config['mapping']['odir'],
        tmp_dir = config['tmpdir'],
        N = "{yid}.%s.{sid}" % config['sambamba_sort']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['sambamba_sort']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['sambamba_sort']['id']),
        j = lambda w: get_resource(w, config, 'sambamba_sort')
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'sambamba_sort')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/sambamba_sort.py"

rule bam_stat:
    input: "{yid}/%s/{sid}.bam" % config['mapping']['odir']
    output: "{yid}/%s/{sid}.tsv" % config['mapping']['odir']
    params:
        N = "{yid}.%s.{sid}" % config['bam_stat']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bam_stat']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bam_stat']['id']),
        j = lambda w: get_resource(w, config, 'bam_stat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bam_stat')['ppn']
    conda: "../envs/work.yml"
    shell: "bam.py stat {input} > {output}"

def merge_bamstats_inputs(w):
    yid = w.yid
    inputs = []
    for sid in config['y'][yid]['SampleID']:
        inputs.append("%s/%s/%s.tsv" % (yid, config['mapping']['odir'], sid))
    return inputs

rule merge_bamstats:
    input:
        lambda w: expand("%s/%s/{sid}.tsv" % (w.yid, config['mapping']['odir']), sid = config['y'][w.yid]['SampleID'])
    output: protected("%s/{yid}/%s" % (config['oid'], config['merge_bamstats']['out']))
    params:
        N = "{yid}.%s" % (config['merge_bamstats']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['merge_bamstats']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['merge_bamstats']['id']),
        j = lambda w: get_resource(w, config, 'merge_bamstats'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'merge_bamstats')['ppn']
    conda: "../envs/work.yml"
    shell: "merge.stats.R --opt bam_stat -o {output} {input}"



