rule wgc1_prepT:
    input:
        fna = lambda w: config['g'][w.genome]['ref'],
    output:
        expand("%s/{{genome}}/11_chroms/{tchrom}.fa" % config['wgc']['dir1'], \
                tchrom = config['g']['B73']['chroms'].split())
    params:
        odir = "%s/{genome}" % config['wgc']['dir1'],
        cdir = "%s/{genome}/11_chroms" % config['wgc']['dir1'],
        N = "%s.{genome}" % config['wgc']['prepT']['id'],
        e = "%s/%s/{genome}.e" % (config['dirj'], config['wgc']['prepT']['id']),
        o = "%s/%s/{genome}.o" % (config['dirj'], config['wgc']['prepT']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepT')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepT')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepT')['mem']
    threads: config['wgc']['prepT']['ppn']
    conda: "../envs/python.yml"
    shell:
        """
        mkdir -p {params.cdir}
        faSplit byname {input.fna} {params.cdir}/
        """

rule wgc1_prepQ:
    input:
        fna = lambda w: config['g'][w.genome]['ref'],
        chrom_size = lambda w: config['g'][w.genome]['chrom_size'],
        chrom_bed = lambda w: config['g'][w.genome]['chrom_bed'],
        gap = lambda w: config['g'][w.genome]['gap'],
    output:
        fnas = expand("%s/{{genome}}/87_qrys/part.{idx}.fna" % \
                config['wgc']['dir1'],
                idx = range(1,config['wgc']['npieces']+1)),
        chain = "%s/{genome}/86.chain" % config['wgc']['dir1'],
    params:
        odir = "%s/{genome}" % config['wgc']['dir1'],
        cdir = "%s/{genome}/11_chroms" % config['wgc']['dir1'],
        npieces = config['wgc']['npieces'],
        N = "%s.{genome}" % config['wgc']['prepQ']['id'],
        e = "%s/%s/{genome}.e" % (config['dirj'], config['wgc']['prepQ']['id']),
        o = "%s/%s/{genome}.o" % (config['dirj'], config['wgc']['prepQ']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepQ')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepQ')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepQ')['mem']
    threads: config['wgc']['prepQ']['ppn']
    conda: "../envs/python.yml"
    shell:
        """
        bed.py filter -min 5000 {input.gap} > {params.odir}/81.qry.gap.bed
        subtractBed -nonamecheck -sorted -a {input.chrom_bed} \
            -b {input.gap} | bed.py filter -min 100 - | \
            bed.py makewindow -w 1000000 -s 995000 - \
            > {params.odir}/85.qry.clean.bed
        bed.py size {params.odir}/85.qry.clean.bed
        bed.py binpacking {params.odir}/85.qry.clean.bed {params.odir}/86.bed {params.odir}/87_qrys --N {params.npieces}
        cut -f1,3 {params.odir}/86.bed > {params.odir}/86.sizes
        chain.py fromBed {params.odir}/86.bed {params.odir}/86.sizes {input.chrom_size} > {output.chain}
        ls {params.odir}/87_qrys/*.bed | parallel "fasta.py extract {input.fna} {{}} > {{.}}.fna"
        rm {params.odir}/87_qrys/*.bed
        """
#        faSplit.R {params.odir}/85.qry.clean.bed \
#            {input.fna} {input.chrom_size} {params.odir}/87_qrys \
#            --chain {params.odir}/86.chain -n {params.npieces}

rule wgc2_align:
    input:
        tgt = "%s/{tgt}/11_chroms/{tchrom}.fa" % config['wgc']['dir1'],
        qry = "%s/{qry}/87_qrys/part.{idx}.fna" % config['wgc']['dir1'],
    output:
        "%s/{qry}_{tgt}/q{idx}.{tchrom}.psl" % config['wgc']['dir2']
    params:
        sam = "%s/{qry}_{tgt}/q{idx}.{tchrom}.sam" % config['wgc']['dir2'],
        lav = "%s/{qry}_{tgt}/q{idx}.{tchrom}.lav" % config['wgc']['dir2'],
        N = lambda w: "%s.%s%s.%s.%s" % (config['wgc']['align']['id'], w.qry[:1], w.tgt[:1], w.idx, w.tchrom),
        e = lambda w: "%s/%s/%s.%s/%s.%s.e" % (config['dirj'], config['wgc']['align']['id'], w.qry, w.tgt, w.idx, w.tchrom),
        o = lambda w: "%s/%s/%s.%s/%s.%s.o" % (config['dirj'], config['wgc']['align']['id'], w.qry, w.tgt, w.idx, w.tchrom),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'align')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'align')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'align')['mem']
    threads: config['wgc']['align']['ppn']
    conda: "../envs/python.yml"
    shell:
        """
        minimap2 -x asm20 -a -t {threads} \
                {input.tgt} {input.qry} | \
                sam.py 2psl - {output}
        """
        #lastz {input.tgt} {input.qry} \
        #        --show=defaults --progress \
        #        --gfextend --gapped --inner=1000 --format=lav > {params.lav}
        #lavToPsl {params.lav} {output}
        #pblat {input.tgt_bit} {input.qry} -threads={threads} \
        #        -ooc={input.tgt_ooc} {output}
        #rm {params.lav}

rule wgc3_merge:
    input:
        expand("%s/{{qry}}_{{tgt}}/q{idx}.{tchrom}.psl" % \
                config['wgc']['dir2'], \
                idx = range(1,config['wgc']['npieces']+1), \
                tchrom = config['g']['B73']['chroms'].split())
    output:
        "%s/{qry}_{tgt}/02.coord.pass.psl" % config['wgc']['dir3']
    params:
        odir = "%s/{qry}_{tgt}" % config['wgc']['dir3'],
        chain = lambda w: "%s/%s/86.chain" %
            (config['wgc']['dir1'], w.qry),
        tbit = lambda w: config['g'][w.tgt]['blat']['x.2bit'],
        qbit = lambda w: config['g'][w.qry]['blat']['x.2bit'],
        tsize = lambda w: config['g'][w.tgt]['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['chrom_size'],
        N = "%s.{qry}.{tgt}" % config['wgc']['merge']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc']['merge']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc']['merge']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'merge')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'merge')['mem']
    threads: config['wgc']['merge']['ppn']
    conda: "../envs/python.yml"
    shell:
        """
        pslCat -nohead {input} > {params.odir}/01.psl
        psl.py coordQ {params.odir}/01.psl {params.qsize} \
                >{params.odir}/02.coord.psl
        mkdir -p {params.odir}
        pslCheck {params.odir}/02.coord.psl \
            -pass={params.odir}/02.coord.pass.psl \
            -fail={params.odir}/02.coord.fail.psl || echo non_success
        """
        #pslSwap {params.odir}/01.psl {params.odir}/02.swap.psl
        #liftOver -pslT {params.odir}/02.swap.psl {params.chain} \
        #        {params.odir}/03.swap.coord.psl {params.odir}/unMapped
        #pslCheck -querySizes={params.qsize} -targetSizes={params.tsize} \
        #        {params.odir}/04.psl \
        #        -pass={params.odir}/05.pass.psl \
        #        -fail={params.odir}/05.fail.psl

rule wgc4_chain:
    input:
        "%s/{qry}_{tgt}/02.coord.pass.psl" % config['wgc']['dir3']
    output:
        "%s/{qry}_{tgt}/15.chain" % config['wgc']['dir3']
    params:
        odir = "%s/{qry}_{tgt}" % config['wgc']['dir3'],
        tbit = lambda w: config['g'][w.tgt]['blat']['x.2bit'],
        qbit = lambda w: config['g'][w.qry]['blat']['x.2bit'],
        tsize = lambda w: config['g'][w.tgt]['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['chrom_size'],
        N = "%s.{qry}.{tgt}" % config['wgc']['chain']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc']['chain']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc']['chain']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'chain')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'chain')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'chain')['mem']
    threads: config['wgc']['chain']['ppn']
    conda: "../envs/python.yml"
    shell:
        """
        axtChain -linearGap=medium -psl {input} \
                {params.tbit} {params.qbit} {params.odir}/10.chain
        chainPreNet {params.odir}/10.chain {params.tsize} {params.qsize} \
                {params.odir}/11.chain
        chainSwap {params.odir}/11.chain {params.odir}/11.q.chain

        chainNet {params.odir}/11.chain {params.tsize} {params.qsize} \
                {params.odir}/13.t.net {params.odir}/13.q.net
        netChainSubset {params.odir}/13.t.net {params.odir}/11.chain stdout | \
                chainSort stdin {params.odir}/13.t.chain
        netChainSubset {params.odir}/13.q.net {params.odir}/11.q.chain stdout | \
                chainSort stdin {params.odir}/13.q.chain
        chainNet {params.odir}/13.q.chain {params.qsize} {params.tsize} \
                /dev/null {params.odir}/15.net
        netChainSubset {params.odir}/15.net {params.odir}/11.chain \
                {params.odir}/15.chain
        """

rule wgc5_post: # python+gatk+R
    input:
        "%s/{qry}_{tgt}/15.chain" % config['wgc']['dir3']
    output:
        "%s/{qry}_{tgt}/10.vnt.bed" % config['dirr'],
        "%s/{qry}_{tgt}/10.{tgt}.vcf.gz" % config['dirr'],
        "%s/{qry}_{tgt}/10.{qry}.vcf.gz" % config['dirr'],
    params:
        odir = "%s/{qry}_{tgt}" % config['dirr'],
        tpre = "%s/{qry}_{tgt}/10.{tgt}" % config['dirr'],
        qpre = "%s/{qry}_{tgt}/10.{qry}" % config['dirr'],
        tfas = lambda w: config['g'][w.tgt]['ref'],
        qfas = lambda w: config['g'][w.qry]['ref'],
        tsize = lambda w: config['g'][w.tgt]['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['chrom_size'],
        tref = lambda w: config['g'][w.tgt]['gatk']['xref'],
        qref = lambda w: config['g'][w.qry]['gatk']['xref'],
        N = "%s.{qry}.{tgt}" % config['wgc']['post']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc']['post']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc']['post']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'post')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'post')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'post')['mem']
    threads: config['wgc']['post']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        chainStitchId {input} {params.odir}/01.stitched.chain
        chainFilter -minGapless=1000 {params.odir}/01.stitched.chain \
            > {params.odir}/02.filtered.chain
        chain.py 2bed {params.odir}/02.filtered.chain > {params.odir}/02.bed

        chainBedVnt.R {params.odir}/02.bed {params.odir}/05.itv.bed
        wgc.py callvnt {params.odir}/05.itv.bed {params.tfas} {params.qfas} \
            --vnt {params.odir}/05.vnt.bed > {params.odir}/05.bed

        chainBedFilter.R {params.odir}/05.bed {params.odir}/10.bed \
            {params.odir}/05.vnt.bed {params.odir}/10.vnt.bed

        chain.py fromBed {params.odir}/10.bed {params.tsize} {params.qsize} \
            > {params.odir}/10.{wildcards.tgt}_{wildcards.qry}.chain
        chainSwap {params.odir}/10.{wildcards.tgt}_{wildcards.qry}.chain \
            {params.odir}/10.{wildcards.qry}_{wildcards.tgt}.chain

        wgc.py bed2vcf --tgt {wildcards.tgt} --qry {wildcards.qry} \
            {params.odir}/10.vnt.bed \
            {params.tpre}.1.vcf {params.qpre}.1.vcf
        sortBed -header -i {params.tpre}.1.vcf >{params.tpre}.1.s.vcf
        sortBed -header -i {params.qpre}.1.vcf >{params.qpre}.1.s.vcf
        gatk UpdateVCFSequenceDictionary -V \
            {params.tpre}.1.s.vcf -R {params.tref} -O {params.tpre}.2.vcf
        gatk UpdateVCFSequenceDictionary -V \
            {params.qpre}.1.s.vcf -R {params.qref} -O {params.qpre}.2.vcf
        bcftools norm -f {params.tfas} -c w -d all \
            {params.tpre}.2.vcf -Ou |\
            bcftools sort -Oz -o {params.tpre}.vcf.gz
        bcftools norm -f {params.qfas} -c w -d all \
            {params.qpre}.2.vcf -Ou |\
            bcftools sort -Oz -o {params.qpre}.vcf.gz
        rm {params.tpre}.[12].*
        rm {params.qpre}.[12].*
        """

rule wgc6_eff_t: # python+snpeff
    input:
        "%s/{qry}_{tgt}/10.{tgt}.vcf.gz" % config['dirr'],
    output:
        "%s/{qry}_{tgt}/15.{tgt}.tsv" % config['dirr'],
    params:
        odir = "%s/{qry}_{tgt}" % config['dirr'],
        tfas = lambda w: config['g'][w.tgt]['ref'],
        qfas = lambda w: config['g'][w.qry]['ref'],
        tsize = lambda w: config['g'][w.tgt]['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['chrom_size'],
        tref = lambda w: config['g'][w.tgt]['gatk']['xref'],
        qref = lambda w: config['g'][w.qry]['gatk']['xref'],
        txcfg = lambda w: config['g'][w.tgt]['snpeff']['xpre'],
        N = "%s.{qry}.{tgt}.t" % config['wgc']['eff']['id'],
        e = "%s/%s/{qry}.{tgt}.t.e" % (config['dirj'], config['wgc']['eff']['id']),
        o = "%s/%s/{qry}.{tgt}.t.o" % (config['dirj'], config['wgc']['eff']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'snpeff')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'snpeff')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'snpeff')['mem']
    threads: config['snpeff']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.txcfg} {wildcards.tgt} -ud 0 \
                {input[0]} | bgzip > {params.odir}/15.{wildcards.tgt}.vcf.gz
        wgc.py parseEff {params.odir}/15.{wildcards.tgt}.vcf.gz > {output[0]}
        """

rule wgc6_eff_q:
    input:
        "%s/{qry}_{tgt}/10.{qry}.vcf.gz" % config['dirr'],
    output:
        "%s/{qry}_{tgt}/15.{qry}.tsv" % config['dirr'],
    params:
        odir = "%s/{qry}_{tgt}" % config['dirr'],
        tfas = lambda w: config['g'][w.tgt]['ref'],
        qfas = lambda w: config['g'][w.qry]['ref'],
        tsize = lambda w: config['g'][w.tgt]['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['chrom_size'],
        tref = lambda w: config['g'][w.tgt]['gatk']['xref'],
        qref = lambda w: config['g'][w.qry]['gatk']['xref'],
        qxcfg = lambda w: config['g'][w.qry]['snpeff']['xpre'],
        N = "%s.{qry}.{tgt}.q" % config['wgc']['eff']['id'],
        e = "%s/%s/{qry}.{tgt}.q.e" % (config['dirj'], config['wgc']['eff']['id']),
        o = "%s/%s/{qry}.{tgt}.q.o" % (config['dirj'], config['wgc']['eff']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'snpeff')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'snpeff')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'snpeff')['mem']
    threads: config['snpeff']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.qxcfg} {wildcards.qry} -ud 0 \
                {input[0]} | bgzip > {params.odir}/15.{wildcards.qry}.vcf.gz
        wgc.py parseEff {params.odir}/15.{wildcards.qry}.vcf.gz > {output[0]}
        """

rule wgc_orthofinder:
    input:
        tgt_faa = lambda w: config['g'][w.tgt]['lfaa'],
        qry_faa = lambda w: config['g'][w.qry]['lfaa'],
    output:
        "%s/{qry}_{tgt}/10.tsv" % config['wgc']['dir4'],
    params:
        odir = "%s/{qry}_{tgt}" % config['wgc']['dir4'],
        odir1 = "%s/{qry}_{tgt}/01_seqs" % config['wgc']['dir4'],
        tfaa = "%s/{qry}_{tgt}/{tgt}.faa" % config['wgc']['dir4'],
        qfaa = "%s/{qry}_{tgt}/{qry}.faa" % config['wgc']['dir4'],
        N = "%s.{qry}.{tgt}" % config['orthofinder']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['orthofinder']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['orthofinder']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'orthofinder')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'orthofinder')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'orthofinder')['mem']
    threads: config['orthofinder']['ppn']
    conda: "../envs/blast.yml"
    shell:
        """
        mkdir -p {params.odir}
        cp -f {input.tgt_faa} {params.tfaa}
        cp -f {input.qry_faa} {params.qfaa}
        orthofinder -f {params.odir} -t {threads} -p {params.odir} -og
        """

