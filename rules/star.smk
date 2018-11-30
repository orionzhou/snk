def star_inputs(w):
    sid = w.sid
    idir = config['star']['idir']
    inputs = dict()
    if config['t'][sid]['paired']:
        inputs['fq1'] = "%s/%s_1.fq.gz" % (idir, sid)
        inputs['fq2'] = "%s/%s_2.fq.gz" % (idir, sid)
    else:
        inputs['fq'] = "%s/%s.fq.gz" % (idir, sid)
    return inputs

def star_input_str(w):
    sid = w.sid
    idir = config['star']['idir']
    input_str = ''
    if config['t'][sid]['paired']:
        fq1 = "%s/%s_1.fq.gz" % (idir, sid)
        fq2 = "%s/%s_2.fq.gz" % (idir, sid)
        input_str = "%s %s" % (fq1, fq2)
    else:
        fq = "%s/%s.fq.gz" % (idir, sid)
        input_str = "%s" % (fq)
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
        unpack(star_inputs)
    output:
        temp("%s/{sid}/Aligned.out.bam" % config['star']['odir1']),
        "%s/{sid}/Log.final.out" % config['star']['odir1']
    log:
        "%s/%s/{sid}.log" % (config['dirl'], config['star']['id'])
    params:
        index = config[config['reference']]["star"],
        input_str = lambda w: star_input_str(w),
        outprefix = lambda w: "%s/%s/" % (config['star']['odir1'], w.sid),
        readcmd = "--readFilesCommand zcat",
        extra = star_extra,
        N = lambda w: "%s.%s" % (config['star']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['star']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['star']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
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

rule sambamba_sort:
    input:
        "%s/{sid}/Aligned.out.bam" % config['star']['odir1']
    output:
        "%s/{sid}.bam" % config['star']['odir2']
    params:
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra']),
        N = lambda w: "%s.%s" % (config['sambamba']['sort']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['sambamba']['sort']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['sambamba']['sort']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['mem']
    threads: config['sambamba']['ppn']
    shell:
        "sambamba sort {params.extra} -t {threads} -o {output} {input}"



