def bismark_inputs(w):
    sid = w.sid
    idir = config['bismark']['idir']
    inputs = dict()
    if config['t'][sid]['paired']:
        inputs['fq1'] = "%s/%s_1.fq.gz" % (idir, sid)
        inputs['fq2'] = "%s/%s_2.fq.gz" % (idir, sid)
    else:
        inputs['fq'] = "%s/%s.fq.gz" % (idir, sid)
    return inputs

def bismark_input_str(w):
    sid = w.sid
    idir = config['bismark']['idir']
    input_str = ''
    if config['t'][sid]['paired']:
        fq1 = "%s/%s_1.fq.gz" % (idir, sid)
        fq2 = "%s/%s_2.fq.gz" % (idir, sid)
        input_str = "-1 %s -2 %s" % (fq1, fq2)
    else:
        fq = "%s/%s.fq.gz" % (idir, sid)
        input_str = "%s" % (fq)
    return input_str

def bismark_extra(w):
    extras = [config["bismark"]["extra"]]
    sm = w.sid
    if 'Genotype' in config['t'][w.sid]:
        sm = config['t'][w.sid]['Genotype']
    extras.append("--rg_tag --rg_id %s --rg_sample %s" % (w.sid, sm))
    return " ".join(extras)

rule bismark:
    input:
        unpack(bismark_inputs)
    output:
        temp("%s/{sid}.sam" % config['bismark']['odir1']),
    params:
        index = config[config['reference']]["bismark"],
        input_str = lambda w: bismark_input_str(w),
        odir = config['bismark']['odir1'],
        extra = bismark_extra,
        N = lambda w: "%s.%s" % (config['bismark']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bismark']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bismark']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bismark')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bismark')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bismark')['mem']
    threads: config['bismark']['ppn']
    shell:
        """
        bismark --parallel {threads} \
        {params.index} {params.input_str} \
        {params.extra} \
        --output_dir {params.odir} --basename {wildcards.sid}
        """

rule sambamba_sort:
    input:
        "%s/{sid}.bam" % config['bismark']['odir1']
    output:
        "%s/{sid}.bam" % config['bismark']['odir2']
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

rule bam_stat:
    input:
        "%s/{sid}.bam" % config['bismark']['odir2']
    output:
        "%s/{sid}.tsv" % config['bismark']['odir2']
    params:
        N = lambda w: "%s.%s" % (config['bam_stat']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bam_stat']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bam_stat']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bam_stat')['mem']
    threads: config['bam_stat']['ppn']
    shell:
        "bam.py stat {input} > {output}"


