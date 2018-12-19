def bismark_extract_extra(wildcards):
    extra = config["bismark_extract"]["extra"]
    return extra

rule sambamba_sort_name:
    input:
        "%s/{sid}.bam" % config['bismark_extract']['idir']
    output:
        "%s/{sid}.rn.bam" % config['bismark_extract']['idir']
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
        "sambamba sort {params.extra} -n -t {threads} -o {output} {input}"

rule bismark_extract:
    input:
        "%s/{sid}.rn.bam" % config['bismark_extract']['idir']
    output:
        "%s/{sid}.rds" % config['bismark_extract']['odir']
    params:
        index = config[config['reference']]["bismark"],
        odir = config['bismark_extract']['odir'],
        cx_report = lambda w: "%s/%s.rn.CX_report.txt" % (config['bismark_extract']['odir'], w.sid),
        parallel = lambda w, resources: int(resources.ppn / 2),
        N = lambda w: "%s.%s" % (config['bismark_extract']['id'], w.sid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['bismark_extract']['id'], w.sid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['bismark_extract']['id'], w.sid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bismark_extract')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bismark_extract')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bismark_extract')['mem']
    threads: config['bismark_extract']['ppn']
    shell:
        """
        bismark_methylation_extractor --parallel {params.parallel} \
        --buffer_size {params.mem}G \
        --genome_folder {params.index} \
        -p --no_overlap \
        --cytosine_report --CX \
        --output {params.odir} \
        {input}
        bsm2gr.R {params.cx_report} {output}
        """
#--bedGraph --zero_based \
#--comprehensive \
# --no-header

