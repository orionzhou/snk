rule fq_uncompress:
    input:
        r0 = ancient("%s/{yid}/{sid}.fq.gz" % config['trimming']['idir']),
        r1 = ancient("%s/{yid}/{sid}_1.fq.gz" % config['trimming']['idir']),
        r2 = ancient("%s/{yid}/{sid}_2.fq.gz" % config['trimming']['idir']),
    output:
        r0 = temp("{yid}/%s/{sid}.fq" % config['trimming']['od12']),
        r1 = temp("{yid}/%s/{sid}_1.fq" % config['trimming']['od12']),
        r2 = temp("{yid}/%s/{sid}_2.fq" % config['trimming']['od12']),
    params:
        N = "{yid}.%s.{sid}" % config['fq_uncompress']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['fq_uncompress']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['fq_uncompress']['id']),
        j = lambda w: get_resource(w, config, 'fq_uncompress'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fq_uncompress')['ppn']
    run:
        if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
            shell("pigz -p {threads} -cd {input.r1} > {output.r1}")
            shell("pigz -p {threads} -cd {input.r2} > {output.r2}")
            shell("touch {output.r0}")
        else:
            shell("pigz -p {threads} -cd {input.r0} > {output.r0}")
            shell("touch {output.r1} {output.r2}")

def fq_extract_range(w):
    yid, sid, part = w.yid, w.sid, w.part
    i = config['y'][yid]['t'][sid]['p'][part]
    part_size = config['y'][yid]['t'][sid]['part_size']

    lines_per_file = part_size * 4
    line_start = (i - 1) * lines_per_file + 1
    line_end = i * lines_per_file
    return "%d,%d" % (line_start, line_end)

rule fq_extract:
    input:
        r0 = ancient("{yid}/%s/{sid}.fq" % config['trimming']['od12']),
        r1 = ancient("{yid}/%s/{sid}_1.fq" % config['trimming']['od12']),
        r2 = ancient("{yid}/%s/{sid}_2.fq" % config['trimming']['od12']),
    output:
        r0 = temp("{yid}/%s/{sid}_{part}.fq" % config['trimming']['od12']),
        r1 = temp("{yid}/%s/{sid}_{part}_1.fq" % config['trimming']['od12']),
        r2 = temp("{yid}/%s/{sid}_{part}_2.fq" % config['trimming']['od12']),
    params:
        linerange = fq_extract_range,
        N = "{yid}.%s.{sid}{part}" % config['fq_extract']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['fq_extract']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['fq_extract']['id']),
        j = lambda w: get_resource(w, config, 'fq_extract'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fq_extract')['ppn']
    run:
        if config['y'][wildcards.yid]['t'][wildcards.sid]['paired']:
            if config['y'][wildcards.yid]['t'][wildcards.sid]['npart'] == 1 and wildcards.part == 'a':
                shell("ln -f %s {output.r1}" % op.abspath(input.r1))
                shell("ln -f %s {output.r2}" % op.abspath(input.r2))
            else:
                shell("sed -n '{params.linerange} p' {input.r1} > {output.r1}")
                shell("sed -n '{params.linerange} p' {input.r2} > {output.r2}")
            shell("touch {output.r0}")
        else:
            if config['y'][wildcards.yid]['t'][wildcards.sid]['npart'] == 1 and wildcards.part == 'a':
                shell("ln -f %s {output.r0}" % op.abspath(input.r0))
            else:
                shell("sed -n '{params.linerange} p' {input.r0} > {output.r0}")
            shell("touch {output.r1} {output.r2}")


def fastp_extra(w):
    yid, sid, part = w.yid, w.sid, w.part
    readtype = config['y'][yid]['readtype']
    assert readtype in config['valid']['readtype'], "invalid readtype: %s" % readtype
    extra = ''
    if readtype == '3rnaseq':
        extra += '-x '
    return extra

rule fastp:
    input:
        r0 = ancient("{yid}/%s/{sid}_{part}.fq" % config['trimming']['od12']),
        r1 = ancient("{yid}/%s/{sid}_{part}_1.fq" % config['trimming']['od12']),
        r2 = ancient("{yid}/%s/{sid}_{part}_2.fq" % config['trimming']['od12'])
    output:
        r0 = "{yid}/%s/{sid}_{part}.fq.gz" % config['trimming']['od14p'],
        r1 = "{yid}/%s/{sid}_{part}_1.fq.gz" % config['trimming']['od14p'],
        r2 = "{yid}/%s/{sid}_{part}_2.fq.gz" % config['trimming']['od14p'],
        json = "{yid}/%s/{sid}_{part}.json" % config['trimming']['od14p'],
        html = "{yid}/%s/{sid}_{part}.html" % config['trimming']['od14p']
    params:
        paired = lambda w: int(config['y'][w.yid]['t'][w.sid]['paired']),
        extra = fastp_extra,
        N = "{yid}.%s.{sid}{part}" % config['fastp']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['fastp']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['fastp']['id']),
        j = lambda w: get_resource(w, config, 'fastp'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'fastp')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/fastp.py"

rule trimmomatic:
    input:
        r0 = "%s/{yid}/{sid}_{part}.fq" % config['trimming']['od12'],
        r1 = "%s/{yid}/{sid}_{part}_1.fq" % config['trimming']['od12'],
        r2 = "%s/{yid}/{sid}_{part}_2.fq" % config['trimming']['od12']
    output:
        r0 = "{yid}/%s/{sid}_{part}.fq.gz" % config['trimming']['od14t'],
        r1 = "{yid}/%s/{sid}_{part}_1.fq.gz" % config['trimming']['od14t'],
        r2 = "{yid}/%s/{sid}_{part}_2.fq.gz" % config['trimming']['od14t'],
        r1u = "{yid}/%s/{sid}_{part}_1.unpaired.fq.gz" % config['trimming']['od14t'],
        r2u = "{yid}/%s/{sid}_{part}_2.unpaired.fq.gz" % config['trimming']['od14t']
    params:
        trimmer = [
            "ILLUMINACLIP:%s:2:30:10:8:no" % config['trimmomatic']['adapter_pe'],
            "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:35"],
        N = "{yid}.%s.{sid}{part}" % config['trimmomatic']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['trimmomatic']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['trimmomatic']['id']),
        j = lambda w: get_resource(w, config, 'trimmomatic'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'trimmomatic')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/trimmomatic.py"

rule bbduk:
    input: "{yid}/%s/{sid}_{part}.fq" % config['trimming']['od12']
    output:
        r0 = "{yid}/%s/{sid}_{part}.fq.gz" % config['trimming']['od14b'],
        r1 = "{yid}/%s/{sid}_{part}_1.fq.gz" % config['trimming']['od14b'],
        r2 = "{yid}/%s/{sid}_{part}_2.fq.gz" % config['trimming']['od14b'],
        json = "{yid}/%s/{sid}_{part}.json" % config['trimming']['od14b'],
    params:
        cmd = config['bbduk']['cmd'],
        extra = "ref=%s %s" %
            (','.join(config['bbduk']['refs']), config['bbduk']['extra']),
        N = "{yid}.%s.{sid}{part}" % config['bbduk']['id'],
        e = "{yid}/%s/%s/{sid}{part}.e" % (config['dirj'], config['bbduk']['id']),
        o = "{yid}/%s/%s/{sid}{part}.o" % (config['dirj'], config['bbduk']['id']),
        j = lambda w: get_resource(w, config, 'bbduk'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'bbduk')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} in={input} out={output.r0} json=t {params.extra} >{output.json} 2>&1
        touch {output.r1} {output.r2}
        """

def trimming_inputs(w):
    yid, sid, part = w.yid, w.sid, w.part
    readtype = config['y'][yid]['readtype']
    assert readtype in config['valid']['readtype'], "invalid readtype: %s" % readtype
    idir = ''
    if readtype == '3rnaseq':
        idir = config['trimming']['od14b']
    elif readtype in ['illumina','solid']:
        idir = config['trimming']['od14p']
    pre = op.abspath("%s/%s/%s" % (yid, idir, sid))
    pre = "%s/%s/%s_%s" % (yid, idir, sid, part)
    return dict(
        r0 = ancient("%s.fq.gz" % pre),
        r1 = ancient("%s_1.fq.gz" % pre),
        r2 = ancient("%s_2.fq.gz" % pre),
        json = ancient("%s.json" % pre)
    )

rule trimming:
    input: unpack(trimming_inputs)
    output:
        r0 = "{yid}/%s/{sid}_{part}.fq.gz" % config['trimming']['odir'],
        r1 = "{yid}/%s/{sid}_{part}_1.fq.gz" % config['trimming']['odir'],
        r2 = "{yid}/%s/{sid}_{part}_2.fq.gz" % config['trimming']['odir'],
        json = "{yid}/%s/{sid}_{part}.json" % config['trimming']['odir'],
    run:
        shell("ln -f %s %s" % (op.abspath(input.r0), output.r0))
        shell("ln -f %s %s" % (op.abspath(input.r1), output.r1))
        shell("ln -f %s %s" % (op.abspath(input.r2), output.r2))
        shell("ln -f %s %s" % (op.abspath(input.json), output.json))

def merge_trimstats_inputs(w):
    yid = w.yid
    inputs = []
    for sid in config['y'][yid]['SampleID']:
        pair_suf = 'pe' if config['y'][yid]['t'][sid]['paired'] else 'se'
        for part in config['y'][yid]['t'][sid]['parts']:
            inputs.append("%s/%s/%s_%s.json" % (yid, config['trimming']['odir'], sid, part))
    return inputs

rule merge_trimstats:
    input: merge_trimstats_inputs
    output:
        protected("%s/{yid}/%s" % (config['oid'], config['trimming']['out']))
    params:
        opt = lambda w: 'bbduk' if config['y'][w.yid]['readtype'] == '3rnaseq' else 'fastp',
        N = "{yid}.%s" % (config['merge_trimstats']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['merge_trimstats']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['merge_trimstats']['id']),
        j = lambda w: get_resource(w, config, 'merge_trimstats'),
    resources: attempt = lambda w, attempt: attempt
    conda: "../envs/work.yml"
    shell: "jsonutil.py {params.opt} {input} > {output}"

