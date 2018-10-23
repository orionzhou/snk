def multiqc_inputs(wildcards):
    inputs = []
    for sid in config['SampleID']:
        paired = config['t'][sid]['paired']
        
        trimmer = "bbduk" if config['readtype'] == '3rnaseq' else "trimmomatic"
        trimmer_suf = 'pe' if paired else 'se'
        inputs.append("%s/%s_%s/%s.log" % (config['dirl'], config[trimmer]['id'], trimmer_suf, sid))
        
        if config['mapper'] == 'hisat2':
            inputs.append("%s/%s.txt" % (config['hisat2']['odir1'], sid))
        elif config['mapper'] == 'star':
            if paired:
                inputs.append("%s/%s_p/Log.final.out" % (config['star']['odir1'], sid))
                inputs.append("%s/%s_u/Log.final.out" % (config['star']['odir1'], sid))
            else:
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

rule merge_featurecounts:
    input:
        expand(["%s/{sid}.txt" % config['merge_featurecounts']['idir']], sid = config['SampleID'])
    params:
        N = lambda w: "%s" % (config['merge_featurecounts']['id']),
        e = lambda w: "%s/%s.e" % (config['dirp'], config['merge_featurecounts']['id']),
        o = lambda w: "%s/%s.o" % (config['dirp'], config['merge_featurecounts']['id']),
    output:
        protected("%s/%s" % (config['dird'], config['merge_featurecounts']['out']))
    shell:
        "merge.featurecounts.R -o {output} {input}"

def bamstat_dir():
    dirb = ''
    if config['mapper'] in ['star','hisat2']:
        dirb = config[config['mapper']]['odir2']
    elif config['mapper'] in ['bwa']:
        dirb = config['cleanbam']['odir2']
    return dirb 

rule merge_bamstats:
    input:
        expand("%s/{sid}.tsv" % bamstat_dir(), sid = config['SampleID'])
    output:
        protected("%s/%s" % (config['dird'], config['merge_bamstats']['out']))
    shell:
        "merge.bamstats.R -o {output} {input}"

