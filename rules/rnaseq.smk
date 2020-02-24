from jcvi.utils.natsort import natsorted

def featurecounts_extra(w):
    yid = w.yid
    extra = "--primary -Q 20 -t exon -g gene_id --byReadGroup"
    extra += " --tmpDir %s" % config['tmpdir']
    extra += " -p"
    if config['y'][yid]['stranded'] == 'reverse':
        extra += " -s 2"
    elif config['y'][yid]['stranded'] == 'yes':
        extra += " -s 1"
    return extra

rule featurecounts:
    input: "{yid}/%s/{sid}.bam" % config['rnaseq']['idir']
    output:
        "{yid}/%s/{sid}.txt" % config['rnaseq']['od31f'],
        "{yid}/%s/{sid}.txt.summary" % config['rnaseq']['od31f']
    params:
        gtf = lambda w: db_index(w, 'gtf'),
        extra = featurecounts_extra,
        N = "{yid}.%s.{sid}" % config['featurecounts']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['featurecounts']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['featurecounts']['id']),
        j = lambda w: get_resource(w, config, 'featurecounts'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'featurecounts')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        featureCounts -T {threads} {params.extra} \
            -a {params.gtf} -o {output[0]} {input}
        """

rule merge_featurecounts:
    input:
        lambda w: expand("%s/%s/{sid}.txt" % (w.yid, config['rnaseq']['od31f']), sid = config['y'][w.yid]['SampleID'])
    output: protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_fcnt']))
    params:
        N = "{yid}.%s" % (config['merge_featurecounts']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['merge_featurecounts']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['merge_featurecounts']['id']),
        j = lambda w: get_resource(w, config, 'merge_featurecounts'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'merge_featurecounts')['ppn']
    conda: "../envs/r.yml"
    shell: "merge.stats.R --opt featurecounts -o {output} {input}"

def mmquant_strand(w):
    yid = w.yid
    if config['y'][yid]['stranded'] == 'reverse':
        return "RF"
    elif config['y'][yid]['stranded'] == 'yes':
        return "FR"
    else:
        return "U"

rule mmquant:
    input: "{yid}/%s/{sid}.bam" % config['rnaseq']['idir']
    output:
        "{yid}/%s/{sid}.tsv" % config['rnaseq']['od31m'],
        "{yid}/%s/{sid}.txt" % config['rnaseq']['od31m']
    params:
        gtf = lambda w: db_index(w, 'gtf'),
#        gtf = '/home/springer/zhoux379/projects/genome/data/B73xMo17/15_intervals/02.chrom.gtf',
        strand = mmquant_strand,
        N = "{yid}.%s.{sid}" % config['mmquant']['id'],
        e = "{yid}/%s/%s/{sid}.e" % (config['dirj'], config['mmquant']['id']),
        o = "{yid}/%s/%s/{sid}.o" % (config['dirj'], config['mmquant']['id']),
        j = lambda w: get_resource(w, config, 'mmquant'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'mmquant')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        mmquant -a {params.gtf} -r {input} \
            -o {output[0]} -O {output[1]} \
            -l 0.3 -t {threads} -s {params.strand}
        """

rule merge_mmquant:
    input:
        lambda w: expand("%s/%s/{sid}.tsv" % (w.yid, config['rnaseq']['od31m']), sid = config['y'][w.yid]['SampleID'])
    output: protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_mmq']))
    params:
        N = "{yid}.%s" % (config['merge_mmquant']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['merge_mmquant']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['merge_mmquant']['id']),
        j = lambda w: get_resource(w, config, 'merge_mmquant'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'merge_mmquant')['ppn']
    conda: "../envs/r.yml"
    shell: "merge.stats.R --opt mmquant -o {output} {input}"

rule rc2cpm:
    input:
        sam = lambda w: config['y'][w.yid]['samplelist'],
        exp = "%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_fcnt']),
        cfg = lambda w: config['g'][config['y'][w.yid]['ref']]["rds"]['xout']
    output: protected("%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_cpm']))
    params:
        N = "{yid}.%s" % (config['rc2cpm']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['rc2cpm']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['rc2cpm']['id']),
        j = lambda w: get_resource(w, config, 'rc2cpm'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'rc2cpm')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        rc2cpm.R {input.sam} {input.exp} {output} \
            --opt featurecounts --yid {wildcards.yid} --config {input.cfg}
        """

rule rc2cpm2:
    input:
        sam = "%s/11_qc/{yid}/%s" % (config['dirh'], config['rnaseq']['meta']),
        exp = "%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_fcnt']),
        cfg = lambda w: config['g'][config['y'][w.yid]['ref']]["rds"]['xout']
    output: protected("%s/11_qc/{yid}/%s" % (config['dirh'], config['rnaseq']['out_cpm']))
    params:
        N = "{yid}.%s" % (config['rc2cpm']['id']),
        e = "{yid}/%s/%s.e" % (config['dirj'], config['rc2cpm']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['rc2cpm']['id']),
        j = lambda w: get_resource(w, config, 'rc2cpm'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'rc2cpm')['ppn']
    conda: "../envs/r.yml"
    shell:
        """
        rc2cpm.R {input.sam} {input.exp} {output} \
            --opt featurecounts --yid {wildcards.yid} --config {input.cfg}
        """

rule ril1_gt:
    input:
        lambda w: expand("%s/%s/{s}.bam" % (w.yid, config['rnaseq']['idir']), s = config['y'][w.yid]['SampleID'])
    output:
        vcf = "{yid}/%s/{rid}.vcf.gz" % config['rnaseq']['od41'],
        tbi = "{yid}/%s/{rid}.vcf.gz.tbi" % config['rnaseq']['od41']
    params:
        ref = lambda w: db_index(w, 'gatk'),
        vnt = lambda w: config['ril_variant'][w.yid],
        reg = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        N = "{yid}.%s.{rid}" % config['ril1_gt']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['ril1_gt']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['ril1_gt']['id']),
        j = lambda w: get_resource(w, config, 'ril1_gt'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ril1_gt')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        bcftools mpileup -f {params.ref} -r {params.reg} -Ou {input} |\
            bcftools call -m -C alleles -T {params.vnt} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

def ril2_inputs(w):
    yid = w.yid
    rids = natsorted(config['g'][config['y'][yid]['ref']]['win56'].keys())
    vcfs = expand("%s/%s/{rid}.vcf.gz" %
        (yid, config['rnaseq']['od41']), rid = rids)
    tbis = expand("%s/%s/{rid}.vcf.gz.tbi" %
        (yid, config['rnaseq']['od41']), rid = rids)
    return {'vcfs':vcfs, 'tbis':tbis}

rule ril2_merge_gt:
    input: unpack(ril2_inputs)
    output: "{yid}/%s" % config['rnaseq']['of42']
    params:
        N = "{yid}.%s" % config['ril2_merge_gt']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['ril2_merge_gt']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['ril2_merge_gt']['id']),
        j = lambda w: get_resource(w, config, 'ril2_merge_gt'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ril2_merge_gt')['ppn']
    conda: "../envs/work.yml"
    shell: """
    bcftools concat -n {input.vcfs} -Oz -o {output}
    bcftools index -t {output}
    """

rule ril3_snpbinner:
    input: "{yid}/%s" % config['rnaseq']['of42']
    output: "{yid}/%s/{rid}.csv" % config['rnaseq']['od43']
    params:
        o1 = "{yid}/%s/{rid}.1.tsv" % config['rnaseq']['od43'],
        o2 = "{yid}/%s/{rid}.2.txt" % config['rnaseq']['od43'],
        min_ratio = 0.0001,
        min_bin_size = 100000,
        reg = lambda w: config['g'][config['y'][w.yid]['ref']]['win11'][w.rid],
        N = "{yid}.%s.{rid}" % config['ril3_snpbinner']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['ril3_snpbinner']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['ril3_snpbinner']['id']),
        j = lambda w: get_resource(w, config, 'ril3_snpbinner'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ril3_snpbinner')['ppn']
    conda: "../envs/ril.yml"
    shell: """
    (bcftools query -l {input} | tr "\\n" "\\t" |\
      sed "s/^/marker\\tposition\\t/; s/\\t$/\\n/" &&
    bcftools query -i 'INFO/DP>=3' -r {params.reg} -f'%CHROM\\t%POS[\\t%GT]\\n' {input} |\
      sed -e 's/\.\/\./\-/g; s/0\/0/a/g; s/1\/1/b/g; s/0\/1/h/g' |\
      awk "{{OFS=\\"\\t\\"}}; {{\$1 = \$1\\"_\\"\$2; print}}") > {params.o1}

    snpbinner crosspoints -i {params.o1} -o {params.o2} -r {params.min_ratio}
    snpbinner bins -i {params.o2} -o {output} -l {params.min_bin_size}
    """

rule ril4_mergebin:
    input:
        lambda w: expand("%s/%s/{rid}.csv" % (w.yid, config['rnaseq']['od43']),
                rid = natsorted(config['g'][config['y'][w.yid]['ref']]['win11'].keys())),
    output: "%s/{yid}/%s" % (config['oid'], config['rnaseq']['out_ril'])
    params:
        N = "{yid}.%s" % config['ril4_mergebin']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['ril4_mergebin']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['ril4_mergebin']['id']),
        j = lambda w: get_resource(w, config, 'ril4_mergebin'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'ril4_mergebin')['ppn']
    conda: "../envs/r.yml"
    shell: "merge.stats.R --opt snpbinner -o {output} {input}"

