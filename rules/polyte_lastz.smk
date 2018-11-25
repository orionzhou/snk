rule lastz:
    input:
        "%s/{genotype}_{tgt}.txt" % config['lastz']['idir']
    output:
        protected("%s/{genotype}_{tgt}.tsv" % config['lastz']['odir'])
    params:
        dummy = "%s/dummy" % config['lastz']['idir'],
        odir = lambda w: "%s/%s_%s" % (config['lastz']['odir'], w.genotype, w.tgt),
        fcmd = lambda w: "%s/%s_%s.sh" % (config['lastz']['odir'], w.genotype, w.tgt),
        N = lambda w: "%s.%s.%s" % (config['lastz']['id'], w.genotype, w.tgt),
        e = lambda w: "%s/%s/%s_%s.e" % (config['dirp'], config['lastz']['id'], w.genotype, w.tgt),
        o = lambda w: "%s/%s/%s_%s.o" % (config['dirp'], config['lastz']['id'], w.genotype, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'lastz')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'lastz')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'lastz')['mem']
    threads: config['lastz']["ppn"]
    shell:
        """
        mkdir -p {params.odir}
        create-lastz-jobs.py {wildcards.genotype} {wildcards.tgt} \
                {input} {params.odir} > {params.fcmd}
        bash {params.fcmd}
        parse-lav.py {params.odir} {input} {params.dummy} > {output}
        """


