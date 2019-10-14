rule ph1a_filter:
    input: lambda w: config['y'][w.yid]['vcf']
    output:
        snp = "{yid}/%s/{rid}.snp.vcf.gz" % config['phylo']['od10'],
        idl = "{yid}/%s/{rid}.idl.vcf.gz" % config['phylo']['od10'],
        snp_stat = "{yid}/%s/{rid}.snp.txt" % config['phylo']['od10'],
        idl_stat = "{yid}/%s/{rid}.idl.txt" % config['phylo']['od10'],
    params:
        min_gq = 30,
        min_prop = 0.2,
        min_maf = 0.05,
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        samples = lambda w: ",".join(config['y'][w.yid]['samples']),
        tmp = config['tmpdir'],
        snp_sites = lambda w, output: output.snp.replace(".vcf.gz", ".sites.vcf.gz"),
        idl_sites = lambda w, output: output.idl.replace(".vcf.gz", ".sites.vcf.gz"),
        snp_bed = lambda w, output: output.snp.replace(".vcf.gz", ".bed"),
        idl_bed = lambda w, output: output.idl.replace(".vcf.gz", ".bed"),
        N = "{yid}.%s.{rid}" % config['ph1a_filter']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['ph1a_filter']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['ph1a_filter']['id']),
        j = lambda w: get_resource(w, config, 'ph1a_filter'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ph1a_filter')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        bcftools view -a -s {params.samples} -r {params.region} {input} -Ou | \
          bcftools view -v snps -m2 -M2 -i 'F_PASS(GQ>={params.min_gq} & GT!="mis") >= {params.min_prop} && N_PASS(GT="AA") > 0 && MAF>={params.min_maf}' \
          -Oz -o {output.snp}
        bcftools view -a -s {params.samples} -r {params.region} {input} -Ou | \
          bcftools view -v indels -m2 -M2 -i 'F_PASS(GQ>={params.min_gq} & GT!="mis") >= {params.min_prop} && N_PASS(GT="AA") > 0 && MAF>={params.min_maf}' \
          -Oz -o {output.idl}
        bcftools index -t {output.snp}
        bcftools index -t {output.idl}
        bcftools stats -s - {output.snp} > {output.snp_stat}
        bcftools stats -s - {output.idl} > {output.idl_stat}
        bcftools view -G {output.snp} -Ou | bcftools annotate -x INFO -Oz -o {params.snp_sites}
        bcftools view -G {output.idl} -Ou | bcftools annotate -x INFO -Oz -o {params.idl_sites}
        bcftools index -t {params.snp_sites}
        bcftools index -t {params.idl_sites}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.snp} -o {params.snp_bed}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.idl} -o {params.idl_bed}
        """

rule ph1b_concat:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.snp.vcf.gz" % (w.yid, config['phylo']['od10']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
        stats = lambda w: expand("%s/%s/{rid}.snp.txt" % (w.yid, config['phylo']['od10']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
    output:
        vcf = "{yid}/%s" % config['phylo']['of11'],
        stat = "{yid}/%s" % config['phylo']['of11b'],
    params:
        tmp = config['tmpdir'],
        N = "{yid}.%s" % config['ph1b_concat']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['ph1b_concat']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['ph1b_concat']['id']),
        j = lambda w: get_resource(w, config, 'ph1b_concat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ph1b_concat')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        bcftools concat -Oz -o {output.vcf} {input.vcfs}
        bcftools index -t {output.vcf}
        merge.bcftools.stats.R --simple {input.stats} {output.stat}
        """

rule ph2_subsample:
    input:
        vcf = "{yid}/%s" % config['phylo']['of11'],
        stat = "{yid}/%s" % config['phylo']['of11b'],
    output: "{yid}/%s" % config['phylo']['of31']
    params:
        nsites = 1000000,
        refn = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        stat = lambda w, output: str(output).replace(".vcf.gz", ".stats.txt"),
        tmp = config['tmpdir'],
        N = "{yid}.%s" % config['ph2_subsample']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['ph2_subsample']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['ph2_subsample']['id']),
        j = lambda w: get_resource(w, config, 'ph2_subsample'),
        mem = lambda w: get_resource(w, config, 'ph2_subsample')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ph2_subsample')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/vcf_subsample.py"

rule ph3_iqtree:
    input: "{yid}/%s" % config['phylo']['of31']
    output: "%s/{yid}/%s" % (config['oid'], config['phylo']['of35'])
    params:
        phy = lambda w, output: str(output).replace(".nwk", ".phy"),
        vphy = lambda w, output: str(output).replace(".nwk", ".varsites.phy"),
        pre = lambda w, output: str(output).replace(".nwk", ""),
        tree = lambda w, output: str(output).replace(".nwk", ".treefile"),
        tmp = config['tmpdir'],
        N = "{yid}.%s" % config['ph3_iqtree']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['ph3_iqtree']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['ph3_iqtree']['id']),
        j = lambda w: get_resource(w, config, 'ph3_iqtree'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ph3_iqtree')['ppn']
    conda: "../envs/work.yml"
    shell:
#iqtree -s {params.phy} -m GTR+I+G -bb 1000 -nt {threads} -pre {params.pre}
#iqtree -s {params.phy} -m MFP+ASC -bb 1000 -nt {threads} -pre {params.pre}
        """
        vcf2phylip.py -i {input} -o {params.phy}
        iqtree -s {params.phy} -m HKY85 -bb 1000 -nt {threads} -pre {params.pre}
        cp {params.tree} {output}
        """


