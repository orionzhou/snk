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
        cx = "%s/{yid}/{mid}.cx.gz" % config['oid'],
    params:
        cx = "%s/{yid}/{mid}.cx" % config['oid'],
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
        """

rule bs4_intersect:
    input: "%s/{yid}/{mid}.cx.gz" % config['oid']
    output: "%s/{yid}/{mid}.bed.gz" % config['oid'],
    params:
        itv = lambda w: config['g'][config['y'][w.yid]['ref']]["annotation"]['bsseq'],
        o1a = "%s/{yid}/{mid}.1.cg.bed" % config['oid'],
        o1b = "%s/{yid}/{mid}.1.chg.bed" % config['oid'],
        o1c = "%s/{yid}/{mid}.1.chh.bed" % config['oid'],
        o2a = "%s/{yid}/{mid}.2.cg.bed" % config['oid'],
        o2b = "%s/{yid}/{mid}.2.chg.bed" % config['oid'],
        o2c = "%s/{yid}/{mid}.2.chh.bed" % config['oid'],
        out = "%s/{yid}/{mid}.bed" % config['oid'],
        N = "{yid}.%s.{mid}" % config['bs4_intersect']['id'],
        e = "{yid}/%s/%s/{mid}.e" % (config['dirj'], config['bs4_intersect']['id']),
        o = "{yid}/%s/%s/{mid}.o" % (config['dirj'], config['bs4_intersect']['id']),
        j = lambda w: get_resource(w, config, 'bs4_intersect'),
        mem = lambda w: get_resource(w, config, 'bs4_intersect')['mem'] - 5
    resources: attempt = lambda w, attempt: attempt
    threads: config['bs4_intersect']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        zcat {input} | bioawk -t \
            '{{if($4+$5 >= 0 && $6=="CG") print $1,$2-1,$2,$3,$4,$5,$6}}' \
            > {params.o1a}
        zcat {input} | bioawk -t \
            '{{if($4+$5 >= 0 && $6=="CHG") print $1,$2-1,$2,$3,$4,$5,$6}}' \
            > {params.o1b}
        zcat {input} | bioawk -t \
            '{{if($4+$5 >= 0 && $6=="CHH") print $1,$2-1,$2,$3,$4,$5,$6}}' \
            > {params.o1c}

        intersectBed -a {params.itv} -b {params.o1a} -wo -sorted | \
            groupBy -i stdin -g 1-6 -c 11,11,12,13 -o count,sum,sum,distinct > {params.o2a}
        intersectBed -a {params.itv} -b {params.o1b} -wo -sorted | \
            groupBy -i stdin -g 1-6 -c 11,11,12,13 -o count,sum,sum,distinct > {params.o2b}
        intersectBed -a {params.itv} -b {params.o1c} -wo -sorted | \
            groupBy -i stdin -g 1-6 -c 11,11,12,13 -o count,sum,sum,distinct > {params.o2c}
        cat {params.o2a} {params.o2b} {params.o2c} | \
            sort -S {params.mem}G -k1,1 -k2,2n -k3,3n > {params.out}
        gzip {params.out}
        rm {params.o1a} {params.o1b} {params.o1c}
        rm {params.o2a} {params.o2b} {params.o2c}
        """


