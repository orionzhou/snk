rule filter:
    input:
        lambda w: config['y'][w.yid]['vcf']
    output:
        snp = "{yid}/%s" % config['phylo']['of21'],
        idl = "{yid}/%s" % config['phylo']['of22'],
    params:
        samples = lambda w: ",".join(config['y'][w.yid]['SampleID']),
        tmp = config['tmpdir'],
        snp_stat = lambda w, output: output.snp.replace(".vcf.gz", ".stats.txt"),
        idl_stat = lambda w, output: output.idl.replace(".vcf.gz", ".stats.txt"),
        snp_sites = lambda w, output: output.snp.replace(".vcf.gz", ".sites.vcf.gz"),
        idl_sites = lambda w, output: output.idl.replace(".vcf.gz", ".sites.vcf.gz"),
        snp_bed = lambda w, output: output.snp.replace(".vcf.gz", ".bed"),
        idl_bed = lambda w, output: output.idl.replace(".vcf.gz", ".bed"),
        N = config['phylo_filter']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['phylo_filter']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['phylo_filter']['id']),
        j = lambda w: get_resource(w, config, 'phylo_filter'),
        mem = lambda w: get_resource(w, config, 'phylo_filter')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'phylo_filter')['ppn']
    conda: "../envs/work.yml"
    shell:
#        bcftools view -m2 -M2 -v indels -i 'MAF>0.05' {input} -Oz -o {output.idl}
        """
        bcftools view -m2 -M2 -v snps -i 'MAF>0.05' -s {params.samples} {input} -Oz -o {output.snp}
        bcftools view -m2 -c 4 -v indels -s {params.samples} {input} -Oz -o {output.idl}
        bcftools index -t {output.snp}
        bcftools index -t {output.idl}
        bcftools stats -s - {output.snp} > {params.snp_stat}
        bcftools stats -s - {output.idl} > {params.idl_stat}
        bcftools view -G {output.snp} -Ou | bcftools annotate -x INFO -Oz -o {params.snp_sites}
        bcftools view -G {output.idl} -Ou | bcftools annotate -x INFO -Oz -o {params.idl_sites}
        bcftools index -t {params.snp_sites}
        bcftools index -t {params.idl_sites}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.snp} -o {params.snp_bed}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.idl} -o {params.idl_bed}
        """

rule sample:
    input: "{yid}/%s" % config['phylo']['of21']
    output: "{yid}/%s" % config['phylo']['of31']
    params:
        fraction = 0.01,
        refn = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        stat = lambda w, output: str(output).replace(".vcf.gz", ".stats.txt"),
        tmp = config['tmpdir'],
        N = "%s.{yid}" % config['phylo_sample']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['phylo_sample']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['phylo_sample']['id']),
        j = lambda w: get_resource(w, config, 'phylo_sample'),
        mem = lambda w: get_resource(w, config, 'phylo_sample')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'phylo_sample')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        gatk --java-options "-Xmx{params.mem}G" SelectVariants \
            --tmp-dir {params.tmp} \
            --use-jdk-deflater --use-jdk-inflater \
            -R {params.refn} \
            -V {input} -O {output} -fraction {params.fraction}
        bcftools stats -s - {output} > {params.stat}
        """

rule iqtree:
    input: "{yid}/%s" % config['phylo']['of31']
    output: "%s/{yid}/%s" % (config['oid'], config['phylo']['of35'])
    params:
        phy = lambda w, output: str(output).replace(".treefile", ".phy"),
        vphy = lambda w, output: str(output).replace(".treefile", ".varsites.phy"),
        pre = lambda w, output: str(output).replace(".treefile", ""),
        tmp = config['tmpdir'],
        N = "{yid}.%s" % config['iqtree']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['iqtree']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['iqtree']['id']),
        j = lambda w: get_resource(w, config, 'iqtree'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'iqtree')['ppn']
    conda: "../envs/work.yml"
    shell:
#iqtree -s {params.phy} -m GTR+I+G -bb 1000 -nt {threads} -pre {params.pre}
#iqtree -s {params.phy} -m MFP+ASC -bb 1000 -nt {threads} -pre {params.pre}
#vcf2phylip.py -i {input} -o {params.phy}
        """
        iqtree -s {params.phy} -m HKY85 -bb 1000 -nt {threads} -pre {params.pre}
        """


