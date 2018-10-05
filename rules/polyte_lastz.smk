rule lastz:
    input:
        "%s/{genotype}_{tgt}.txt" % config['lastz']['idir']
    output:
        protected("%s/{genotype}_{tgt}.tsv" % config['lastz']['odir'])
    params:
        dummy = "%s/dummy" % config['lastz']['idir'],
        odir = lambda wildcards: "%s/%s_%s" % (config['lastz']['odir'], wildcards.genotype, wildcards.tgt),
        fcmd = lambda wildcards: "%s/%s_%s.sh" % (config['lastz']['odir'], wildcards.genotype, wildcards.tgt),
    shell:
        """
        mkdir -p {params.odir}
        create-lastz-jobs.py {wildcards.genotype} {wildcards.tgt} \
                {input} {params.odir} > {params.fcmd}
        bash {params.fcmd}
        parse-lav.py {params.odir} {input} {params.dummy} > {output}
        """


