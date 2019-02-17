def mapping_inputs(w):
    yid, sid = w.yid, w.sid
    idir = config['mapping']['idir']
    inputs = dict()
    if config['y'][yid]['t'][sid]['paired']:
        inputs['fq1'] = "%s/%s/%s_1.fq.gz" % (yid, idir, sid)
        inputs['fq2'] = "%s/%s/%s_2.fq.gz" % (yid, idir, sid)
    else:
        inputs['fq'] = "%s/%s/%s.fq.gz" % (yid, idir, sid)
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

def star_extra(w):
    extras = [config["star"]["extra"]]
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
    input:
        unpack(mapping_inputs)
    output:
        temp("{yid}/%s/{sid}/Aligned.out.bam" % config['star']['odir']),
        "{yid}/%s/{sid}/Log.final.out" % config['star']['odir']
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['star']['id'])
    params:
        index = lambda w: config[config['y'][w.yid]['reference']]["star"],
        input_str = lambda w, input: mapping_input_str(w, input, 'star'),
        outprefix = lambda w: "%s/%s/%s/" % (w.yid, config['star']['odir'], w.sid),
        readcmd = "--readFilesCommand zcat",
        extra = star_extra,
        N = "{yid}.%s.{sid}" % config['star']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['star']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['star']['id']),
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'star')['q'],
        ppn = lambda w, attempt:  get_resource(config, attempt, 'star')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'star')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'star')['mem']
    threads: config['star']['ppn']
    shell:
        """
        STAR {params.extra} --runThreadN {threads} \
        --genomeDir {params.index} \
        --readFilesIn {params.input_str} {params.readcmd} \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.outprefix} \
        --outStd Log \
        >{log} 2>&1
        """

def hisat2_extra(w):
    extras = [config["hisat2"]["extra"]]
    extras.append("--rg-id %s --rg SM:%s" % (w.sid, w.sid))
    return " ".join(extras)

rule hisat2:
    input:
        unpack(mapping_inputs)
    output:
        temp("{yid}/%s/{sid}.bam" % config['hisat2']['odir']),
        "{yid}/%s/{sid}.txt" % config['hisat2']['odir']
    params:
        index = lambda w: config[config['y'][w.yid]['reference']]["hisat2"],
        input_str = lambda w, input: mapping_input_str(w, input, 'hisat2'),
        extra = hisat2_extra,
        N = "{yid}.%s.{sid}" % config['hisat2']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['hisat2']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['hisat2']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'hisat2')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'hisat2')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'hisat2')['mem']
    threads: config['hisat2']['ppn']
    shell:
        """
        hisat2 {params.extra} --threads {threads} \
            -x {params.index} {params.input_str} \
            --summary-file {output[1]} \
            | samtools view -Sbh -o {output[0]} -
        """

def bwa_extra(w):
    extras = [config["bwa"]["extra"]]
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
    output:
        temp("{yid}/%s/{sid}.sam" % config['bwa']['odir'])
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['bwa']['id'])
    params:
        index = lambda w: config[config['y'][w.yid]['reference']]["bwa"],
        input_str = lambda w, input: mapping_input_str(w, input, 'bwa'),
        extra = bwa_extra,
        N = "{yid}.%s.{sid}" % config['bwa']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['bwa']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['bwa']['id']),
    resources:
        q = lambda w, attempt:  get_resource(config, attempt, 'bwa')['q'],
        ppn = lambda w, attempt: get_resource(config, attempt, 'bwa')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'bwa')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'bwa')['mem']
    threads: config['bwa']['ppn']
    shell:
        """
        bwa mem -t {threads} {params.index} {params.extra} {params.input_str} \
                >{output} 2>>{log}
        """

def bismark_extra(w):
    extras = [config["bismark"]["extra"]]
    yid, sid = w.yid, w.sid
    sm = sid
    if 'Genotype' in config['y'][yid]['t'][sid]:
        sm = config['y'][yid]['t'][w.sid]['Genotype']
    extras.append("--temp_dir %s" % config['tmpdir'])
    extras.append("--rg_tag --rg_id %s --rg_sample %s" % (sid, sm))
    return " ".join(extras)

rule bismark:
    input:
        unpack(mapping_inputs)
    output:
        bam = "{yid}/%s/{sid}.bam" % config['bismark']['odir'],
        report = "{yid}/%s/{sid}.txt" % config['bismark']['odir'],
    params:
        index = lambda w: config[config['y'][w.yid]['reference']]["bismark"],
        input_str = lambda w, input: mapping_input_str(w, input, 'bismark'),
        odir = config['bismark']['odir'],
        parallel = lambda w, resources: int(resources.ppn / 2),
        extra = bismark_extra,
        N = "{yid}.%s.{sid}" % config['bismark']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['bismark']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['bismark']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bismark')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bismark')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bismark')['mem']
    threads: config['bismark']['ppn']
    shell:
#        --basename {wildcards.sid}
        """
        bismark --parallel {params.parallel} \
            {params.index} {params.input_str} \
            {params.extra} --output_dir {params.odir}
        mv {params.odir}/{wildcards.sid}*bismark*.bam {output.bam}
        mv {params.odir}/{wildcards.sid}*bismark*_report.txt {output.report}
        """

def sambamba_sort_input(w):
    yid, sid = w.yid, w.sid
    mapper = config['y'][yid]['mapper']
    assert mapper in config['valid']['mapper'], "invalid mapper: %s" % mapper
    idir = "%s/%s" % (yid, config[mapper]['odir'])
    fi = ''
    if mapper == 'star':
        fi = "%s/%s/Aligned.out.bam" % (idir, sid)
    elif mapper in ['hisat2','bismark']:
        fi = "%s/%s.bam" % (idir, sid)
    elif mapper == 'bwa':
        fi = "%s/%s.sam" % (idir, sid)
    return fi

rule sambamba_sort:
    input: sambamba_sort_input
    output:
        "{yid}/%s/{sid}.bam" % config['mapping']['odir'],
        "{yid}/%s/{sid}.bam.bai" % config['mapping']['odir']
    params:
        mapper = lambda w: config['y'][w.yid]['mapper'],
        tmp_bam = lambda w, output: "%s.tmp.bam" % output[0],
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra']),
        N = "{yid}.%s.{sid}" % config['sambamba']['sort']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['sambamba']['sort']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['sambamba']['sort']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['mem']
    threads: config['sambamba']['ppn']
    run:
        if params.mapper == 'bwa':
            shell("""
            sambamba view -S -f bam -t {threads} {input} -o {output[0]}
            sambamba sort {params.extra} -t {threads} -o {output[0]} {params.tmp_bam}
            rm {params.tmp_bam}
            """)
        else:
            shell("sambamba sort {params.extra} -t {threads} -o {output[0]} {input}")

rule bam_stat:
    input:
        "{yid}/%s/{sid}.bam" % config['mapping']['odir']
    output:
        "{yid}/%s/{sid}.tsv" % config['mapping']['odir']
    params:
        N = "{yid}.%s.{sid}" % config['bam_stat']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirp'], config['bam_stat']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirp'], config['bam_stat']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['mem']
    threads: config['bam_stat']['ppn']
    shell:
        "bam.py stat {input} > {output}"

def merge_bamstats_inputs(w):
    yid = w.yid
    inputs = []
    for sid in config['y'][yid]['SampleID']:
        inputs.append("%s/%s/%s.tsv" % (yid, config['mapping']['odir'], sid))
    return inputs

rule merge_bamstats:
    input:
        lambda w: expand("%s/%s/{sid}.tsv" % (w.yid, config['mapping']['odir']), sid = config['y'][w.yid]['SampleID'])
    output:
        protected("{yid}/%s/%s" % (config['dird'], config['merge_bamstats']['out']))
    shell:
        "merge.bamstats.R -o {output} {input}"



