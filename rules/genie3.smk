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
        N = lambda w: "%s.%s" % (config['genie3']['id'], w.nid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['genie3']['id'], w.nid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['genie3']['id'], w.nid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'genie3')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'genie3')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'genie3')['mem']
    threads: config['genie3']['ppn']
    shell:
        "genie3.py -p {threads} {params.extra} {input} {output}"


