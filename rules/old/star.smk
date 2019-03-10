def star_inputs(w):
    yid, sid = w.yid, w.sid
    idir = config['star']['idir']
    inputs = dict()
    if config['y'][yid]['t'][sid]['paired']:
        inputs['fq1'] = "%s/%s/%s_1.fq.gz" % (yid, idir, sid)
        inputs['fq2'] = "%s/%s/%s_2.fq.gz" % (yid, idir, sid)
    else:
        inputs['fq'] = "%s/%s/%s.fq.gz" % (yid, idir, sid)
    return inputs

def star_input_str(w, input):
    yid, sid = w.yid, w.sid
    input_str = ''
    if config['y'][yid]['t'][sid]['paired']:
        input_str = "%s %s" % (input['fq1'], input['fq2'])
    else:
        input_str = "%s" % (input['fq'])
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
        temp("{yid}/%s/{sid}/Aligned.out.bam" % config['star']['odir1']),
        "{yid}/%s/{sid}/Log.final.out" % config['star']['odir1']
    log:
        "{yid}/%s/%s/{sid}.log" % (config['dirl'], config['star']['id'])
    params:
        index = config[config['reference']]["star"],
        input_str = lambda w, input: star_input_str(w, input),
        outprefix = lambda w: "%s/%s/%s/" % (w.yid, config['star']['odir1'], w.sid),
        readcmd = "--readFilesCommand zcat",
        extra = star_extra,
        N = lambda w: "%s.%s.%s" % (w.yid, config['star']['id'], w.sid),
        e = lambda w: "%s/%s/%s/%s.e" % (w.yid, config['dirp'], config['star']['id'], w.sid),
        o = lambda w: "%s/%s/%s/%s.o" % (w.yid, config['dirp'], config['star']['id'], w.sid),
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



