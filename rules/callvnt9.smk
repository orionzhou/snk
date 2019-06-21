
rule gatk_genomics_dbimport:
    input: unpack(genomics_dbimport_inputs)
    output: "{yid}/%s/{rid}.cp" % config['callvnt2']['od10']
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win127'][w.rid],
        gvcfs = lambda w, input: ["-V %s" % x for x in input.vcfs],
        odir = "{yid}/%s/{rid}" % config['callvnt2']['od10'],
        batch_size = 100,
        extra = gatk_extra(picard = False, jdk = True),
        N = "{yid}.%s.{rid}" % config['gatk_genomics_dbimport']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['gatk_genomics_dbimport']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['gatk_genomics_dbimport']['id']),
        j = lambda w: get_resource(w, config, 'gatk_genomics_dbimport'),
        mem = lambda w: get_resource(w, config, 'gatk_genomics_dbimport')['mem'] - 5,
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'gatk_genomics_dbimport')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        rm -rf {params.odir}
        {params.cmd} --java-options "-Xmx{params.mem}G" GenomicsDBImport \
            {params.extra} -L {params.region} \
            --genomicsdb-workspace-path {params.odir} \
            --batch-size {params.batch_size} --reader-threads {threads} \
            {params.gvcfs}
        touch {output}
        """

