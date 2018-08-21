def genie3_extra(wildcards):
    nid = wildcards.nid
    extra = ''
    if config['t'][nid]['study'] == 'kremling2018':
        extra = '--cpm'
    return extra

rule genie3:
    input:
        "%s/{nid}.rda" % config['genie3']['idir']
    output:
        "%s/{nid}.rda" % config['genie3']['odir']
    params:
        genie3_extra
    threads:
        config['genie3']['threads']
    shell:
        "genie3.R -p {threads} {params} {input} {output}"


