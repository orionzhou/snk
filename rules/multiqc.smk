def multiqc_inputs(wildcards):
    inputs = []
    for sid in config['SampleID']:
        paired = config['t'][sid]['paired']
        suf = 'pe' if paired else 'se'
        
        trimmer = "bbduk" if config['readtype'] == '3rnaseq' else "trimmomatic"
        inputs.append("%s/%s_%s/%s.log" % (config['dirl'], trimmer, suf, sid))
        
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
        "%s/multiqc.log" % config['dirl']
    params:
        outdir = lambda w, output: op.dirname(output[0]),
        outfile = lambda w, output: op.basename(output[0]),
        extra = '',
        N = lambda w: "multiqc",
        ppn = config['multiqc']['ppn'],
        walltime = config['multiqc']['walltime'],
        mem = config['multiqc']['mem']
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
    output:
        protected("%s/%s" % (config['dird'], config['merge_featurecounts']['out']))
    shell:
        "merge.featurecounts.R -o {output} {input}"

rule merge_bamstats:
    input:
        expand(["%s/{sid}.tsv" % config['merge_bamstats']['idir']], sid = config['SampleID'])
    output:
        protected("%s/%s" % (config['dird'], config['merge_bamstats']['out']))
    shell:
        "merge.bamstats.R -o {output} {input}"

