def genie3_extra(wildcards):
    nid = wildcards.nid
    extra = ''
    #if config['t'][nid]['study'] == 'kremling2018':
    #    extra = '--cpm'
    return extra

rule genie3:
    input:
        "%s/{nid}.pkl" % config['genie3']['idir']
    output:
        "%s/{nid}.pkl" % config['genie3']['odir']
    params:
        extra = genie3_extra,
        ppn = config['genie3']['ppn'],
        walltime = config['genie3']['walltime'],
        mem = config['genie3']['mem']
    threads: config['genie3']['ppn']
    shell:
        "genie3.py -p {threads} {params.extra} {input} {output}"


