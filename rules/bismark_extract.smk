def bismark_extract_extra(w):
    extra = config["bismark_extract"]["extra"]
    return extra

rule sambamba_sort_name:
    input:
        "{yid}/%s/{sid}.bam" % config['bismark_extract']['idir']
    output:
        "{yid}/%s/{sid}.rn.bam" % config['bismark_extract']['idir']
    params:
        extra = "--tmpdir=%s %s" % (config['tmpdir'], config['sambamba']['sort']['extra']),
        N = "{yid}.%s.{sid}" % (config['sambamba']['sort']['id'], w.sid),
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['sambamba']['sort']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['sambamba']['sort']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'sambamba', 'sort')['mem']
    threads: config['sambamba']['ppn']
    conda: "../envs/job.yml"
    shell:
        "sambamba sort {params.extra} -n -t {threads} -o {output} {input}"

rule bismark_extract:
    input:
        "{yid}/%s/{sid}.rn.bam" % config['bismark_extract']['idir']
    output:
        "{yid}/%s/{sid}.rds" % config['bismark_extract']['odir']
    params:
        index = config[config['reference']]["bismark"],
        odir = config['bismark_extract']['odir'],
        cx_report = lambda w: "%s/%s.rn.CX_report.txt" % (config['bismark_extract']['odir'], w.sid),
        parallel = lambda w, resources: int(resources.ppn / 2),
        N = "{yid}.%s.{sid}" % config['bismark_extract']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['bismark_extract']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['bismark_extract']['id'])
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bismark_extract')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bismark_extract')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bismark_extract')['mem']
    threads: config['bismark_extract']['ppn']
    conda: "../envs/job.yml"
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

