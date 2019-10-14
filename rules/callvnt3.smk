rule vt1a_filter:
    input: lambda w: config['y'][w.yid]['vcf']
    output:
        snp = "{yid}/%s/{rid}.snp.vcf.gz" % config['callvnt3']['od10'],
        idl = "{yid}/%s/{rid}.idl.vcf.gz" % config['callvnt3']['od10'],
        snp_stat = "{yid}/%s/{rid}.snp.txt" % config['callvnt3']['od10'],
        idl_stat = "{yid}/%s/{rid}.idl.txt" % config['callvnt3']['od10'],
    params:
        min_gq = 30,
        min_prop = 0,
        min_maf = 0,
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        samples = lambda w: ",".join(config['y'][w.yid]['samples']),
        tmp = config['tmpdir'],
        snp_sites = lambda w, output: output.snp.replace(".vcf.gz", ".sites.vcf.gz"),
        idl_sites = lambda w, output: output.idl.replace(".vcf.gz", ".sites.vcf.gz"),
        snp_bed = lambda w, output: output.snp.replace(".vcf.gz", ".bed"),
        idl_bed = lambda w, output: output.idl.replace(".vcf.gz", ".bed"),
        N = "{yid}.%s.{rid}" % config['vt1a_filter']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['vt1a_filter']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['vt1a_filter']['id']),
        j = lambda w: get_resource(w, config, 'vt1a_filter'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'vt1a_filter')['ppn']
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

rule vt1b_concat:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.snp.vcf.gz" % (w.yid, config['callvnt3']['od10']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
        stats = lambda w: expand("%s/%s/{rid}.snp.txt" % (w.yid, config['callvnt3']['od10']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
    output:
        vcf = "{yid}/%s" % config['callvnt3']['of11'],
        stat = "{yid}/%s" % config['callvnt3']['of11b'],
    params:
        tmp = config['tmpdir'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        N = "{yid}.%s" % config['vt1b_concat']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['vt1b_concat']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['vt1b_concat']['id']),
        j = lambda w: get_resource(w, config, 'vt1b_concat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'vt1b_concat')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        bcftools concat -Ou {input.vcfs} |\
          bcftools norm -f {params.ref} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        merge.bcftools.stats.R --simple {input.stats} {output.stat}
        """

rule vt3_ase:
    input:
        vcf = "{yid}/%s" % config['callvnt3']['of11'],
    output:
        vcf = "{yid}/%s/{gt}.vcf.gz" % config['callvnt3']['od21'],
        tbi = "{yid}/%s/{gt}.vcf.gz.tbi" % config['callvnt3']['od21'],
        stat = "{yid}/%s/{gt}.txt" % config['callvnt3']['od21'],
        vcf2 = "{yid}/%s/{gt}.vcf.gz" % config['callvnt3']['od22'],
        bcf = "{yid}/%s/{gt}.bcf" % config['callvnt3']['od22'],
        csi = "{yid}/%s/{gt}.bcf.csi" % config['callvnt3']['od22'],
    params:
        min_gq = 30,
        tmp = config['tmpdir'],
        vcf2 = "{yid}/%s/{gt}.vcf" % config['callvnt3']['od22'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        N = "{yid}.%s.{gt}" % config['vt3_ase']['id'],
        e = "{yid}/%s/%s/{gt}.e" % (config['dirj'], config['vt3_ase']['id']),
        o = "{yid}/%s/%s/{gt}.o" % (config['dirj'], config['vt3_ase']['id']),
        j = lambda w: get_resource(w, config, 'vt3_ase'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'vt3_ase')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/build_ase_vnt.py"



