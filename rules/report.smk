def multiqc_inputs(w):
    inputs = []
    for sid in config['SampleID']:
        paired = config['t'][sid]['paired']
#        trimmer = "bbduk" if config['readtype'] == '3rnaseq' else "trimmomatic"
#        trimmer_suf = 'pe' if paired else 'se'
#        inputs.append("%s/%s_%s/%s.log" % (config['dirl'], config[trimmer]['id'], trimmer_suf, sid))
        if config['mapper'] == 'hisat2':
            inputs.append("%s/%s.txt" % (config['hisat2']['odir1'], sid))
        elif config['mapper'] == 'star':
            inputs.append("%s/%s/Log.final.out" % (config['star']['odir1'], sid))
        #elif config['mapper'] == 'bwa':
        #    inputs.append("%s/%s.txt" % (config['bwa']['odir2'], sid))
        if config['mapper'] in ['star','hisat2']:
            inputs.append("%s/%s.txt.summary" % (config['featurecounts']['odir'], sid))
    return inputs

rule multiqc:
    input:
        multiqc_inputs
    output:
        "%s/%s" % (config['dird'], config['multiqc']['out']),
    log:
        "%s/%s.log" % (config['dirl'], config['multiqc']['id'])
    params:
        outdir = lambda w, output: op.dirname(output[0]),
        outfile = lambda w, output: op.basename(output[0]),
        extra = '',
        N = lambda w: "%s" % (config['multiqc']['id']),
        e = lambda w: "%s/%s.e" % (config['dirp'], config['multiqc']['id']),
        o = lambda w: "%s/%s.o" % (config['dirp'], config['multiqc']['id']),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'multiqc')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'multiqc')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'multiqc')['mem']
    threads: config['multiqc']['ppn']
    shell:
        "multiqc {params.extra} --force "
        "-o {params.outdir} "
        "-n {params.outfile} "
        "{input} "
        ">{log} 2>&1"

def merge_trimstats_inputs(w):
    inputs = []
    for sid in config['SampleID']:
        pair_suf = 'pe' if config['t'][sid]['paired'] else 'se'
        inputs.append("%s/%s.%s.json" % (config['fastp']['odir'], sid, pair_suf))
    return inputs

rule merge_trimstats:
    input: merge_trimstats_inputs
    output:
        "%s/%s" % (config['dird'], config['merge_trimstats']['out'])
    shell:
        "jsonutil.py fastp {input} > {output}"

rule merge_featurecounts:
    input:
        expand(["%s/{sid}.txt" % config['merge_featurecounts']['idir']], sid = config['SampleID'])
    output:
        "%s/%s" % (config['dird'], config['merge_featurecounts']['out'])
    params:
        N = lambda w: "%s" % (config['merge_featurecounts']['id']),
        e = lambda w: "%s/%s.e" % (config['dirp'], config['merge_featurecounts']['id']),
        o = lambda w: "%s/%s.o" % (config['dirp'], config['merge_featurecounts']['id']),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'merge_featurecounts')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'merge_featurecounts')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'merge_featurecounts')['mem']
    threads: config['merge_featurecounts']['ppn']
    shell:
        "merge.featurecounts.R -o {output} {input}"

def bamstat_dir():
    dirb = ''
    if config['mapper'] in ['star','hisat2']:
        dirb = config[config['mapper']]['odir2']
    elif config['mapper'] in ['bwa']:
        dirb = config['callvnt']['idir']
    elif config['mapper'] in ['bismark']:
        dirb = config['bismark_extract']['idir']
    return dirb

rule bam_stat:
    input:
        "%s/{sid}.bam" % bamstat_dir()
    output:
        "%s/{sid}.tsv" % bamstat_dir()
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

rule sambamba_flagstat:
    input:
        "%s/{sid}.bam" % bamstat_dir()
    output:
        "%s/{sid}.flagstat.tsv" % bamstat_dir()
    params:
        extra = config['sambamba']['flagstat']['extra'],
        N = lambda w: "%s.%s" % (config['sambamba']['flagstat']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['sambamba']['flagstat']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['sambamba']['flagstat']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'flagstat')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'flagstat')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'flagstat')['mem']
    threads: config['sambamba']['ppn']
    shell:
        "sambamba flagstat {params.extra} -t {threads} {input} > {output}"

rule merge_bamstats:
    input:
        expand("%s/{sid}.tsv" % bamstat_dir(), sid = config['SampleID'])
    output:
        "%s/%s" % (config['dird'], config['merge_bamstats']['out'])
    shell:
        "merge.bamstats.R -o {output} {input}"

