def grn1_input(w):
    nid = w.nid
    mid = config['t'][nid]['mid']
    subid = config['t'][nid]['subid']
    diri = "/home/springer/zhoux379/projects/rnaseq/data/raw"
    return "%s/%s/cpm.rds" % (diri, mid)

rule grn1_get_exp_mat:
    input: grn1_input
    output:
        raw = "%s/{nid}.tsv" % config['grn']['od11'],
        filt = "%s/{nid}.tsv" % config['grn']['od12'],
    params:
        subid = lambda w: config['t'][w.nid]['subid'],
        pct_sam_on = 0.05,
        min_var_p = 0.25,
        N = "%s.{nid}" % config['grn1_get_exp_mat']['id'],
        e = "%s/%s/{nid}.e" % (config['dirj'], config['grn1_get_exp_mat']['id']),
        o = "%s/%s/{nid}.o" % (config['dirj'], config['grn1_get_exp_mat']['id']),
        j = lambda w: get_resource(w, config, 'grn1_get_exp_mat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn1_get_exp_mat')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        cpm.filter.R {input} {output.raw} \
            --subid={params.subid} --pct_sam_on 0 --num_sam_on 0
        cpm.filter.R {input} {output.filt} \
            --subid={params.subid} --pct_sam_on {params.pct_sam_on} --min_var_p {params.min_var_p}
        """

def grn2_extra(w):
    nid = w.nid
    extra = ''
#    extra = '--tree_method ET'
#    extra = '--tree_method RF'
#    extra += ' --max_depth 10'
    #if config['t'][nid]['study'] == 'kremling2018':
    #    extra = '--cpm'
    return extra

rule grn2_make:
    input:
        em = "%s/{nid}.tsv" % config['grn']['od12'],
        tf = ancient("%s/%s" % (config['dirh'], config['grn']['tf'])),
    output: temp("%s/{gopt}.{nid}.pkl" % config['grn']['od14'])
    params:
        extra = grn2_extra,
        N = "%s.{gopt}.{nid}" % config['grn2_make']['id'],
        e = "%s/%s/{gopt}.{nid}.e" % (config['dirj'], config['grn2_make']['id']),
        o = "%s/%s/{gopt}.{nid}.o" % (config['dirj'], config['grn2_make']['id']),
        j = lambda w: get_resource(w, config, 'grn2_make'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn2_make')['ppn']
    conda: "../envs/work.yml"
    shell: """
    grn.py --opt {wildcards.gopt} -p {threads} \
        {params.extra} {input.em} {input.tf} {output}"""

rule grn3_convert:
    input:
        em = "%s/{nid}.tsv" % config['grn']['od12'],
        pkl = "%s/{gopt}.{nid}.pkl" % config['grn']['od14']
    output: protected("%s/{gopt}.{nid}.rds" % config['grn']['od15'])
    params:
        N = "%s.{gopt}.{nid}" % config['grn3_convert']['id'],
        e = "%s/%s/{gopt}.{nid}.e" % (config['dirj'], config['grn3_convert']['id']),
        o = "%s/%s/{gopt}.{nid}.o" % (config['dirj'], config['grn3_convert']['id']),
        j = lambda w: get_resource(w, config, 'grn3_convert'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn3_convert')['ppn']
    conda: "../envs/work.yml"
    shell: "grn.pkl2rds.R {input.em} {input.pkl} {output}"

rule grn4_merge:
    input:
        expand("%s/{{gopt}}.{nid}.rds" % config['grn']['od15'], nid = config['nid'])
    output: "%s/{gopt}.rds" % config['oid'],
    params:
        N = "%s.{gopt}" % config['grn4_merge']['id'],
        e = "%s/%s/{gopt}.e" % (config['dirj'], config['grn4_merge']['id']),
        o = "%s/%s/{gopt}.o" % (config['dirj'], config['grn4_merge']['id']),
        j = lambda w: get_resource(w, config, 'grn4_merge'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn4_merge')['ppn']
    conda: "../envs/work.yml"
    shell: "grn.merge.R {output} --gopt {wildcards.gopt}"

rule grn3_cv:
    input:
        f_cfg = "%s/%s" % (config['dirh'], config['grn']['cfg']),
        a_fi_filt = "%s/{nid}.tsv" % config['grn']['od12'],
    output:
        "%s/{nid}.tsv" % config['grn']['od15']
    params:
        nid = "{nid}",
        a_mdir = "%s/{nid}" % config['grn']['od14'],
        dir_filt = config['grn']['od12'],
        dir_raw = config['grn']['od11'],
        dir_model = config['grn']['od14'],
        N = "%s.{nid}" % config['grn3_cv']['id'],
        e = "%s/%s/{nid}.e" % (config['dirj'], config['grn3_cv']['id']),
        o = "%s/%s/{nid}.o" % (config['dirj'], config['grn3_cv']['id']),
        j = lambda w: get_resource(w, config, 'grn3_cv'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn3_cv')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        grn.py eval_model {params.nid} {output} -p {threads} --f_cfg {input.f_cfg} \
            --dir_filt {params.dir_filt} --dir_raw {params.dir_raw} --dir_model {params.dir_model}
        """

rule grn4_merge_cv:
    input:
        expand("%s/{nid}.tsv" % config['grn']['od15'], nid=config['nid']),
    output: "%s/01.meval.rds" % config['oid'],
    params:
        N = "%s" % config['grn4_merge_cv']['id'],
        e = "%s/%s.e" % (config['dirj'], config['grn4_merge_cv']['id']),
        o = "%s/%s.o" % (config['dirj'], config['grn4_merge_cv']['id']),
        j = lambda w: get_resource(w, config, 'grn4_merge_cv'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn4_merge_cv')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        grn.meval.merge.R -o {output} {input}
        """

rule grn7_eval:
    input: "%s/{gopt}.{nid}.rds" % config['grn']['od15']
    output: "%s/{eopt}.{gopt}.{nid}.rds" % config['grn']['od17']
    params:
        N = "%s.{eopt}.{gopt}.{nid}" % config['grn7_eval']['id'],
        e = "%s/%s/{eopt}/{gopt}.{nid}.e" % (config['dirj'], config['grn7_eval']['id']),
        o = "%s/%s/{eopt}/{gopt}.{nid}.o" % (config['dirj'], config['grn7_eval']['id']),
        j = lambda w: get_resource(w, config, 'grn7_eval'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn7_eval')['ppn']
    conda: "../envs/work.yml"
    shell: "grn.eval.R {input} {output} --opt {wildcards.eopt} --permut 1000"

rule grn8_merge_eval:
    input:
        expand("%s/{{eopt}}.{{gopt}}.{nid}.rds" % config['grn']['od17'], nid = config['nid'])
    output: "%s/{gopt}.{eopt}.rds" % config['oid'],
    params:
        N = "%s.{gopt}.{eopt}" % config['grn8_merge_eval']['id'],
        e = "%s/%s/{gopt}.{eopt}.e" % (config['dirj'], config['grn8_merge_eval']['id']),
        o = "%s/%s/{gopt}.{eopt}.o" % (config['dirj'], config['grn8_merge_eval']['id']),
        j = lambda w: get_resource(w, config, 'grn8_merge_eval'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'grn8_merge_eval')['ppn']
    conda: "../envs/work.yml"
    shell: "grn.eval.merge.R {output} --gopt {wildcards.gopt} --eopt {wildcards.eopt}"

