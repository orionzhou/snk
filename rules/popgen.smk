rule pg1_vcf2fas:
    input:
        ref = config['popgen']['ref'],
        vcf = config['popgen']['vcf'],
        cdsloc = config['popgen']['cdsloc']
    output:
        fas = "%s/{sid}.fas" % config['popgen']['od05'],
        fai = "%s/{sid}.fas.fai" % config['popgen']['od05'],
    params:
        N = "%s.{sid}" % config['pg1_vcf2fas']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['pg1_vcf2fas']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['pg1_vcf2fas']['id']),
        j = lambda w: get_resource(w, config, 'pg1_vcf2fas'),
    resources:
        attempt = lambda w, attempt: attempt,
        load = lambda w: get_resource(w, config, 'pg1_vcf2fas')['load'],
    threads: lambda w: get_resource(w, config, 'pg1_vcf2fas')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        samtools faidx {input.ref} -r {input.cdsloc} |\
            bcftools consensus -M N -s {wildcards.sid} {input.vcf} -o {output.fas}
        samtools faidx {output.fas}
        """

rule pg2_cds2gene:
    input:
        gene = config['popgen']['gene'],
        cds = "%s/{sid}.fas" % config['popgen']['od05'],
    output:
        fas = "%s/{sid}.fas" % config['popgen']['od06'],
        fai = "%s/{sid}.fas.fai" % config['popgen']['od06'],
    params:
        pre = "{sid}#",
        N = "%s.{sid}" % config['pg2_cds2gene']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['pg2_cds2gene']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['pg2_cds2gene']['id']),
        j = lambda w: get_resource(w, config, 'pg2_cds2gene'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'pg2_cds2gene')['ppn']
    conda: "../envs/work.yml"
    shell:
        """fasta.py cds2gene --prefix {params.pre} {input.cds} {input.gene} {output.fas}
        samtools faidx {output.fas}"""

rule pg3_merge:
    input: expand("%s/{sid}.fas" % config['popgen']['od06'], sid =config['popgen']['sids'])
    output:
        fas = "%s" % config['popgen']['of10'],
        fai = "%s.fai" % config['popgen']['of10'],
    params:
        N = "%s" % config['pg3_merge']['id'],
        e = "%s/%s.e" % (config['dirj'], config['pg3_merge']['id']),
        o = "%s/%s.o" % (config['dirj'], config['pg3_merge']['id']),
        j = lambda w: get_resource(w, config, 'pg3_merge'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'pg3_merge')['ppn']
    conda: "../envs/work.yml"
    shell:
        """cat {input} > {output.fas}
        samtools faidx {output.fas}"""

rule pg4_stat:
    input:
        sample = config['popgen']['sample_list'],
        gene = config['popgen']['gene_id'],
        fas = "%s" % config['popgen']['of10'],
        fai = "%s.fai" % config['popgen']['of10']
    output: "%s/{batch}.tsv" % config['popgen']['od12']
    params:
        batch_size = config['popgen']['batch_size'],
        max_missing = 0.5,
        N = "%s.{batch}" % config['pg4_stat']['id'],
        e = "%s/%s/{batch}.e" % (config['dirj'], config['pg4_stat']['id']),
        o = "%s/%s/{batch}.o" % (config['dirj'], config['pg4_stat']['id']),
        j = lambda w: get_resource(w, config, 'pg4_stat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'pg4_stat')['ppn']
    conda: "../envs/egglib.yml"
    shell: '''
        popgen.py stats {input.fas} {input.gene} {input.sample} {output} \
        --batch {wildcards.batch} --batch_size {params.batch_size} \
        --max_missing {params.max_missing}'''

rule pg5_merge_stats:
    input: expand("%s/{batch}.tsv" % config['popgen']['od12'], batch=config['popgen']['batches'])
    output: config['popgen']['of13']
    params:
        N = "%s" % config['pg5_merge_stats']['id'],
        e = "%s/%s.e" % (config['dirj'], config['pg5_merge_stats']['id']),
        o = "%s/%s.o" % (config['dirj'], config['pg5_merge_stats']['id']),
        j = lambda w: get_resource(w, config, 'pg5_merge_stats'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'pg5_merge_stats')['ppn']
    shell: "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} > {output}"

