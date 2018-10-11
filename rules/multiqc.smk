rule multiqc:
    input:
        multiqc_inputs
    output:
        "%s/multiqc.html" % config['dird']
    log:
        "%s/multiqc.log" % config['dirl']
    params:
        outdir = lambda wildcards, output: op.dirname(output[0]),
        outfile = lambda wildcards, output: op.basename(output[0]),
        extra = '',
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

