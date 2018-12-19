def genie3_extra(wildcards):
    nid = wildcards.nid
    extra = ''
    #if config['t'][nid]['study'] == 'kremling2018':
    #    extra = '--cpm'
    return extra

rule genie3:
    input:
        "%s/{nid}.pkl" % config['grn']['genie3']['idir']
    output:
        temp("%s/{nid}.pkl" % config['grn']['genie3']['odir'])
    params:
        extra = genie3_extra,
        N = lambda w: "%s.%s" % (config['grn']['genie3']['id'], w.nid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['grn']['genie3']['id'], w.nid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['grn']['genie3']['id'], w.nid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','genie3')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','genie3')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','genie3')['mem']
    threads: config['grn']['genie3']['ppn']
    shell:
        "genie3.py -p {threads} {params.extra} {input} {output}"

rule pkl2rda:
    input:
        "%s/{nid}.pkl" % config['grn']['genie3']['idir'],
        "%s/{nid}.pkl" % config['grn']['pkl2rda']['idir']
    output:
        "%s/{nid}.rda" % config['grn']['pkl2rda']['odir']
    params:
        N = lambda w: "%s.%s" % (config['grn']['pkl2rda']['id'], w.nid),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['grn']['pkl2rda']['id'], w.nid),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['grn']['pkl2rda']['id'], w.nid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rda')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rda')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rda')['mem']
    threads: config['grn']['pkl2rda']['ppn']
    shell:
        "genie3_pkl2rda.R {input[0]} {input[1]} {output}"

rule eval:
    input:
        "%s/{nid}.rda" % config['grn']['eval']['idir']
    output:
        "%s/{nid}_{evtype}.rds" % config['grn']['eval']['odir'],
    params:
        N = lambda w: "%s.%s.%s" % (config['grn']['eval']['id'], w.nid, w.evtype),
        e = lambda w: "%s/%s/%s/%s.e" % (config['dirp'], config['grn']['eval']['id'], w.evtype, w.nid),
        o = lambda w: "%s/%s/%s/%s.o" % (config['dirp'], config['grn']['eval']['id'], w.evtype, w.nid),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval')['mem']
    threads: config['grn']['eval']['ppn']
    shell:
        "grn.eval.R {input} {output} --opt {wildcards.evtype}"

rule eval_merge:
    input:
        expand("%s/{nid}_{{evtype}}.rds" % config['grn']['eval_merge']['idir'], nid = config['nid'])
    output:
        "%s/01.{evtype}.rds" % config['grn']['eval_merge']['odir'],
    params:
        N = lambda w: "%s.%s" % (config['grn']['eval_merge']['id'], w.evtype),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['grn']['eval_merge']['id'], w.evtype),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['grn']['eval_merge']['id'], w.evtype),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['mem']
    threads: config['grn']['eval_merge']['ppn']
    shell:
        "grn.eval.merge.R {output} --opt {wildcards.evtype}"

