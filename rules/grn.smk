def get_exp_mat_input(w):
    nid = w.nid
    mid = config['t'][nid]['mid']
    subid = config['t'][nid]['subid']
    diri = "/home/springer/zhoux379/projects/rnaseq/data/08_raw_output"
    return "%s/%s/cpm.rds" % (diri, mid)

rule get_exp_mat:
    input: get_exp_mat_input
    output:
        "%s/{nid}.tsv" % config['grn']['get_exp_mat']['odir'],
    params:
        subid = lambda w: config['t'][w.nid]['subid'],
        N = "%s.{nid}" % config['grn']['get_exp_mat']['id'],
        e = "%s/%s/{nid}.e" % (config['dirp'], config['grn']['get_exp_mat']['id']),
        o = "%s/%s/{nid}.o" % (config['dirp'], config['grn']['get_exp_mat']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','get_exp_mat')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','get_exp_mat')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','get_exp_mat')['mem']
    threads: config['grn']['get_exp_mat']['ppn']
    shell:
        """
        conda activate r
        cpm.filter.R {input} {output} --subid={params.subid}
        """

def genie3_extra(wildcards):
    nid = wildcards.nid
    extra = ''
    #if config['t'][nid]['study'] == 'kremling2018':
    #    extra = '--cpm'
    return extra

rule genie3:
    input:
        em = "%s/{nid}.tsv" % config['grn']['get_exp_mat']['odir'],
        tf = config['grn']['genie3']['tf']
    output:
        temp("%s/{nid}.pkl" % config['grn']['genie3']['odir'])
    params:
        extra = genie3_extra,
        N = "%s.{nid}" % config['grn']['genie3']['id'],
        e = "%s/%s/{nid}.e" % (config['dirp'], config['grn']['genie3']['id']),
        o = "%s/%s/{nid}.o" % (config['dirp'], config['grn']['genie3']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','genie3')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','genie3')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','genie3')['mem']
    threads: config['grn']['genie3']['ppn']
    shell:
        """
        set +u
        source activate python
        grn.py genie3 -p {threads} {params.extra} {input.em} {input.tf} {output}
        """

rule pkl2rda:
    input:
        em = "%s/{nid}.tsv" % config['grn']['get_exp_mat']['odir'],
        net = "%s/{nid}.pkl" % config['grn']['genie3']['odir']
    output:
        "%s/{nid}.rda" % config['grn']['genie3']['odir']
    params:
        N = "%s.{nid}" % config['grn']['pkl2rda']['id'],
        e = "%s/%s/{nid}.e" % (config['dirp'], config['grn']['pkl2rda']['id']),
        o = "%s/%s/{nid}.o" % (config['dirp'], config['grn']['pkl2rda']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rda')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rda')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rda')['mem']
    threads: config['grn']['pkl2rda']['ppn']
    shell:
        """
        conda activate r
        genie3_pkl2rda.R {input.em} {input.net} {output}
        """

rule eval:
    input:
        "%s/{nid}.rda" % config['grn']['genie3']['odir']
    output:
        "%s/{nid}_{evtype}.rds" % config['grn']['eval']['odir'],
    params:
        N = "%s.{nid}.{evtype}" % config['grn']['eval']['id'],
        e = "%s/%s/{evtype}/{nid}.e" % (config['dirp'], config['grn']['eval']['id']),
        o = "%s/%s/{evtype}/{nid}.o" % (config['dirp'], config['grn']['eval']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval')['mem']
    threads: config['grn']['eval']['ppn']
    shell:
        """
        source activate r
        grn.eval.R {input} {output} --opt {wildcards.evtype}
        """

rule eval_merge:
    input:
        expand("%s/{nid}_{{evtype}}.rds" % config['grn']['eval']['odir'], nid = config['nid'])
    output:
        "%s/01.{evtype}.rds" % config['grn']['eval_merge']['odir'],
    params:
        N = "%s.{evtype}" % config['grn']['eval_merge']['id'],
        e = "%s/%s/{evtype}.e" % (config['dirp'], config['grn']['eval_merge']['id']),
        o = "%s/%s/{evtype}.o" % (config['dirp'], config['grn']['eval_merge']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['mem']
    threads: config['grn']['eval_merge']['ppn']
    shell:
        """
        conda activate r
        grn.eval.merge.R {output} --opt {wildcards.evtype}
        """

