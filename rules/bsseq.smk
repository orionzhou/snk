def bs1_inputs(w):
    yid, mid = w.yid, w.mid
    sids = config['y'][yid]['m'][mid]
    return expand("%s/%s/{sid}.bam" % (yid, config['bsseq']['idir']), sid=sids)

rule bs1_sbb:
    input: bs1_inputs
    output: "{yid}/%s/{mid}.bam" % config['bsseq']['od31']
    params:
        sort = 'byname',
        odir = "{yid}/%s" % config['bsseq']['od31'],
        tmp_dir = config['tmpdir'],
        N = "{yid}.%s.{mid}" % config['sambamba_sort']['id'],
        e = "{yid}/%s/%s/{mid}.e" % (config['dirj'], config['sambamba_sort']['id']),
        o = "{yid}/%s/%s/{mid}.o" % (config['dirj'], config['sambamba_sort']['id']),
        j = lambda w: get_resource(w, config, 'sambamba_sort'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'sambamba_sort')['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/sambamba_sort.py"

def bs2_extra(w):
    yid, mid = w.yid, w.mid
    sids = config['y'][yid]['m'][mid]
    paireds = set([config['y'][yid]['t'][sid]['paired'] for sid in sids])
    if len(paireds) > 1:
        print("mid[ %s ] has both paired and single reads: %s" % (mid, sids))
        return '-p'
#        sys.exit(1)
    elif paireds.pop():
        return '-p'
    else:
        return '-s'

rule bs2_extract:
    input: "{yid}/%s/{mid}.bam" % config['bsseq']['od31']
    output: "{yid}/%s/{mid}.CX_report.txt" % config['bsseq']['od33']
    params:
        index = lambda w: config['g'][config['y'][w.yid]['ref']]["bismark"]['xpre'],
        extra = bs2_extra,
        odir = "{yid}/%s" % config['bsseq']['od33'],
        p = lambda w: int(get_resource(w, config, 'bs2_extract')['ppn'] / 2),
        N = "{yid}.%s.{mid}" % config['bs2_extract']['id'],
        e = "{yid}/%s/%s/{mid}.e" % (config['dirj'], config['bs2_extract']['id']),
        o = "{yid}/%s/%s/{mid}.o" % (config['dirj'], config['bs2_extract']['id']),
        j = lambda w: get_resource(w, config, 'bs2_extract'),
        mem = lambda w: get_resource(w, config, 'bs2_extract')['mem'] - 2
    resources: attempt = lambda w, attempt: attempt
    threads: config['bs2_extract']['ppn']
    conda: "../envs/bismark.yml"
    shell: #--bedGraph --zero_based \ #--comprehensive \ # --no-header
        #--buffer_size {params.mem}G
        """
        bismark_methylation_extractor --parallel {params.p} \
            --ample_memory \
            --genome_folder {params.index} \
            {params.extra} --no_overlap \
            --cytosine_report --CX \
            --output {params.odir} \
            {input}
        rm {params.odir}/*_{wildcards.mid}.txt
        """

rule bs3_convert:
    input: "{yid}/%s/{mid}.CX_report.txt" % config['bsseq']['od33']
    output:
        cx = "%s/{yid}/{mid}.cx.gz" % config['bsseq']['odir'],
        bed1 = "%s/{yid}/{mid}.cg.bed.gz" % config['bsseq']['odir'],
        bed2 = "%s/{yid}/{mid}.chg.bed.gz" % config['bsseq']['odir'],
        bed3 = "%s/{yid}/{mid}.chh.bed.gz" % config['bsseq']['odir'],
        out = "%s/{yid}/{mid}.bed.gz" % config['oid'],
    params:
        itv = lambda w: config['g'][config['y'][w.yid]['ref']]["annotation"]['bsseq'],
        cx = "%s/{yid}/{mid}.cx" % config['bsseq']['odir'],
        bed1 = "%s/{yid}/{mid}.cg.bed" % config['bsseq']['odir'],
        bed2 = "%s/{yid}/{mid}.chg.bed" % config['bsseq']['odir'],
        bed3 = "%s/{yid}/{mid}.chh.bed" % config['bsseq']['odir'],
        out1 = "%s/{yid}/{mid}.cg.bed" % config['oid'],
        out2 = "%s/{yid}/{mid}.chg.bed" % config['oid'],
        out3 = "%s/{yid}/{mid}.chh.bed" % config['oid'],
        out = "%s/{yid}/{mid}.bed" % config['oid'],
        N = "{yid}.%s.{mid}" % config['bs3_convert']['id'],
        e = "{yid}/%s/%s/{mid}.e" % (config['dirj'], config['bs3_convert']['id']),
        o = "{yid}/%s/%s/{mid}.o" % (config['dirj'], config['bs3_convert']['id']),
        j = lambda w: get_resource(w, config, 'bs3_convert'),
        mem = lambda w: get_resource(w, config, 'bs3_convert')['mem'] - 5
    resources: attempt = lambda w, attempt: attempt
    threads: config['bs3_convert']['ppn']
    conda: "../envs/work.yml"
    shell:
#cgmaptools convert bismark2cgmap -i {params.cx} -o {params.cgmap}
#gzip {params.cgmap}
        """
        sort -S {params.mem}G -k1,1 -k2,2n {input} > {params.cx}
        gzip {params.cx}

        zcat {output.cx} | bioawk -t \
            '{{if($4+$5 > 0 && $6=="CG") print $1,$2-1,$2,$3,$4,$5,$6}}' \
            > {params.bed1}
        zcat {output.cx} | bioawk -t \
            '{{if($4+$5 > 0 && $6=="CHG") print $1,$2-1,$2,$3,$4,$5,$6}}' \
            > {params.bed2}
        zcat {output.cx} | bioawk -t \
            '{{if($4+$5 > 0 && $6=="CHH") print $1,$2-1,$2,$3,$4,$5,$6}}' \
            > {params.bed3}

        intersectBed -a {params.itv} -b {params.bed1} -wo -sorted | \
            groupBy -i stdin -g 1-6 -c 11,11,12,13 -o count,sum,sum,distinct > {params.out1}
        intersectBed -a {params.itv} -b {params.bed2} -wo -sorted | \
            groupBy -i stdin -g 1-6 -c 11,11,12,13 -o count,sum,sum,distinct > {params.out2}
        intersectBed -a {params.itv} -b {params.bed3} -wo -sorted | \
            groupBy -i stdin -g 1-6 -c 11,11,12,13 -o count,sum,sum,distinct > {params.out3}
        cat {params.out1} {params.out2} {params.out3} | \
            sort -S {params.mem}G -k1,1 -k2,2n -k3,3n > {params.out}
        gzip {params.bed1}
        gzip {params.bed2}
        gzip {params.bed3}
        gzip {params.out}
        rm {params.out1} {params.out2} {params.out3}
        """


