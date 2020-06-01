rule vt1a_filter:
    input: lambda w: config['y'][w.yid]['vcf']
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['callvnt3']['od10'],
        stat = "{yid}/%s/{rid}.txt" % config['callvnt3']['od10'],
    params:
        min_gq = 30,
        min_prop = 0,
        min_maf = 0,
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        samples = lambda w: ",".join(config['y'][w.yid]['samples']),
        tmp = config['tmpdir'],
        sites = lambda w, output: output.vcf.replace(".vcf.gz", ".sites.vcf.gz"),
        bed = lambda w, output: output.vcf.replace(".vcf.gz", ".bed"),
        N = "{yid}.%s.{rid}" % config['vt1a_filter']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['vt1a_filter']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['vt1a_filter']['id']),
        j = lambda w: get_resource(w, config, 'vt1a_filter'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'vt1a_filter')['ppn']
    conda: "../envs/callvnt.yml"
    shell:
        """
        bcftools view -a -s {params.samples} -r {params.region} {input} -Ou | \
          bcftools view -m2 -i 'ALT!="*" && F_PASS(GQ>={params.min_gq} & GT!="mis") >= {params.min_prop} && N_PASS(GT="AA") > 0 && MAF>={params.min_maf}' \
          -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        bcftools stats -s - {output.vcf} > {output.stat}
        bcftools view -G {output.vcf} -Ou | bcftools annotate -x INFO -Oz -o {params.sites}
        bcftools index -t {params.sites}
        bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' {output.vcf} -o {params.bed}
        """

rule vt1b_concat:
    input:
        vcfs = lambda w: expand("%s/%s/{rid}.vcf.gz" % (w.yid, config['callvnt3']['od10']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win56'].keys())),
        stats = lambda w: expand("%s/%s/{rid}.txt" % (w.yid, config['callvnt3']['od10']),
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
    conda: "../envs/callvnt.yml"
    shell:
        """
        bcftools concat -Ou {input.vcfs} |\
          bcftools norm -f {params.ref} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        merge.bcftools.stats.R --simple {input.stats} {output.stat}
        """

rule vt3_ase:
    input: vcf = "{yid}/%s" % config['callvnt3']['of11'],
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
    conda: "../envs/callvnt.yml"
    script: "../scripts/build_ase_vnt.py"



