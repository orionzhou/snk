def get_exp_mat_input(w):
    nid = w.nid
    mid = config['t'][nid]['mid']
    subid = config['t'][nid]['subid']
    diri = "/home/springer/zhoux379/projects/rnaseq/data/08_raw_output"
    return "%s/%s/cpm.rds" % (diri, mid)

rule get_exp_mat:
    input: get_exp_mat_input
    output:
        raw = "%s/{nid}.tsv" % config['grn']['od11'],
        filt = "%s/{nid}.tsv" % config['grn']['od12'],
    params:
        subid = lambda w: config['t'][w.nid]['subid'],
        pct_sam_on = 0.05,
        min_var_p = 0.25,
        N = "%s.{nid}" % config['grn']['get_exp_mat']['id'],
        e = "%s/%s/{nid}.e" % (config['dirj'], config['grn']['get_exp_mat']['id']),
        o = "%s/%s/{nid}.o" % (config['dirj'], config['grn']['get_exp_mat']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','get_exp_mat')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','get_exp_mat')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','get_exp_mat')['mem']
    threads: config['grn']['get_exp_mat']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        cpm.filter.R {input} {output.raw} \
            --subid={params.subid} --pct_sam_on 0 --num_sam_on 0
        cpm.filter.R {input} {output.filt} \
            --subid={params.subid} --pct_sam_on {params.pct_sam_on} --min_var_p {params.min_var_p}
        """

def make_grn_extra(w):
    nid = w.nid
    extra = '--tree_method ET'
    extra = '--tree_method RF'
    extra += ' --max_depth 10'
    #if config['t'][nid]['study'] == 'kremling2018':
    #    extra = '--cpm'
    return extra

rule make_grn:
    input:
        em = "%s/{nid}.tsv" % config['grn']['od12'],
        tf = "%s/%s" % (config['dird'], config['grn']['tf'])
    output:
        pkl = temp("%s/{nid}.pkl" % config['grn']['od14']),
        odir = directory("%s/{nid}" % config['grn']['od14']),
    params:
        extra = make_grn_extra,
        odir = "%s/{nid}" % config['grn']['od14'],
        N = "%s.{nid}" % config['grn']['make_grn']['id'],
        e = "%s/%s/{nid}.e" % (config['dirj'], config['grn']['make_grn']['id']),
        o = "%s/%s/{nid}.o" % (config['dirj'], config['grn']['make_grn']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','make_grn')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','make_grn')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','make_grn')['mem']
    threads: config['grn']['make_grn']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        grn.py genie3 -p {threads} {params.extra} \
            {input.em} {input.tf} {output.pkl} {output.odir}
        """

rule pkl2rds:
    input:
        em = "%s/{nid}.tsv" % config['grn']['od12'],
        net = "%s/{nid}.pkl" % config['grn']['od14']
    output: "%s/{nid}.rds" % config['grn']['od14']
    params:
        N = "%s.{nid}" % config['grn']['pkl2rds']['id'],
        e = "%s/%s/{nid}.e" % (config['dirj'], config['grn']['pkl2rds']['id']),
        o = "%s/%s/{nid}.o" % (config['dirj'], config['grn']['pkl2rds']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rds')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rds')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','pkl2rds')['mem']
    threads: config['grn']['pkl2rds']['ppn']
    conda: "../envs/work.yml"
    shell:
        "grn.pkl2rds.R {input.em} {input.net} {output}"

rule eval_model:
    input:
        f_cfg = op.join(config['dird'], config['grn']['cfg']),
        a_fi_filt = "%s/{nid}.tsv" % config['grn']['od12'],
        a_mdir = "%s/{nid}" % config['grn']['od14'],
    output:
        "%s/{nid}.tsv" % config['grn']['od15']
    params:
        nid = "{nid}",
        dir_filt = config['grn']['od12'],
        dir_raw = config['grn']['od11'],
        dir_model = config['grn']['od14'],
        N = "%s.{nid}" % config['grn']['eval_model']['id'],
        e = "%s/%s/{nid}.e" % (config['dirj'], config['grn']['eval_model']['id']),
        o = "%s/%s/{nid}.o" % (config['dirj'], config['grn']['eval_model']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval_model')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval_model')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval_model')['mem']
    threads: config['grn']['eval_model']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        grn.py eval_model {params.nid} {output} -p {threads} --f_cfg {input.f_cfg} \
            --dir_filt {params.dir_filt} --dir_raw {params.dir_raw} --dir_model {params.dir_model}
        """

rule merge_eval_model:
    input:
        expand("%s/{nid}.tsv" % config['grn']['od15'], nid=config['nid']),
    output: "%s/01.meval.rds" % config['dirr'],
    params:
        N = "%s" % config['grn']['merge_eval_model']['id'],
        e = "%s/%s.e" % (config['dirj'], config['grn']['merge_eval_model']['id']),
        o = "%s/%s.o" % (config['dirj'], config['grn']['merge_eval_model']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','merge_eval_model')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','merge_eval_model')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','merge_eval_model')['mem']
    threads: config['grn']['merge_eval_model']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        grn.meval.merge.R -o {output} {input}
        """

rule eval:
    input: "%s/{nid}.rds" % config['grn']['od14']
    output: "%s/{nid}_{evtype}.rds" % config['grn']['od17']
    params:
        N = "%s.{nid}.{evtype}" % config['grn']['eval']['id'],
        e = "%s/%s/{evtype}/{nid}.e" % (config['dirj'], config['grn']['eval']['id']),
        o = "%s/%s/{evtype}/{nid}.o" % (config['dirj'], config['grn']['eval']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval')['mem']
    threads: config['grn']['eval']['ppn']
    conda: "../envs/work.yml"
    shell:
        "grn.eval.R {input} {output} --opt {wildcards.evtype} --permut 500"

rule eval_merge:
    input:
        expand("%s/{nid}_{{evtype}}.rds" % config['grn']['od17'], nid = config['nid'])
    output: "%s/01.{evtype}.rds" % config['dirr'],
    params:
        N = "%s.{evtype}" % config['grn']['eval_merge']['id'],
        e = "%s/%s/{evtype}.e" % (config['dirj'], config['grn']['eval_merge']['id']),
        o = "%s/%s/{evtype}.o" % (config['dirj'], config['grn']['eval_merge']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'grn','eval_merge')['mem']
    threads: config['grn']['eval_merge']['ppn']
    conda: "../envs/work.yml"
    shell:
        "grn.eval.merge.R {output} --opt {wildcards.evtype}"

