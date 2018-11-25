rule wgc1_prepare:
    input:
        lambda w: config[w.genome]['ref']
    output:
        fna = "%s/{genome}/10_genome.fna" % config['wgc']['dir1'],
        size = "%s/{genome}/15_intervals/01.chrom.sizes" % config['wgc']['dir1'],
    params:
        odir = lambda w: "%s/%s" % (config['wgc']['dir1'], w.genome),
        cdir = lambda w: "%s/%s/11_chroms" % (config['wgc']['dir1'], w.genome),
        extra = '',
        N = lambda w: "%s.%s" % (config['wgc']['prepare']['id'], w.genome),
        e = lambda w: "%s/%s/%s.e" % (config['dirp'], config['wgc']['prepare']['id'], w.genome),
        o = lambda w: "%s/%s/%s.o" % (config['dirp'], config['wgc']['prepare']['id'], w.genome),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepare')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepare')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'prepare')['mem']
    threads: config['wgc']['prepare']['ppn']
    shell:
        """
        mkdir -p {params.odir}
        """

rule wgc1_break_tgt:
    input:
        lambda w: config[w.genome]['ref'],
        #"%s/{genome}/10_genome.fna" % config['wgc']['dir1']
    output:
        expand("%s/{{genome}}/11_chroms/{tchrom}.fa" % config['wgc']['dir1'], \
                tchrom = config['B73']['chroms'].split())
    params:
        odir = lambda w: "%s/%s" % (config['wgc']['dir1'], w.genome),
        cdir = lambda w: "%s/%s/11_chroms" % (config['wgc']['dir1'], w.genome),
    shell:
        """
        mkdir -p {params.cdir}
        faSplit byname {input} {params.cdir}/
        """

rule wgc1_break_qry:
    input:
        lambda w: config[w.genome]['ref'],
        #"%s/{genome}/10_genome.fna" % config['wgc']['dir1']
    output:
        fnas = expand("%s/{{genome}}/87_qrys/part.{idx}.fna" % \
                config['wgc']['dir1'], 
                idx = range(1,config['wgc']['npieces']+1)),
        chain = "%s/{genome}/86.chain" % config['wgc']['dir1'],
    params:
        odir = lambda w: "%s/%s" % (config['wgc']['dir1'], w.genome),
        cdir = lambda w: "%s/%s/11_chroms" % (config['wgc']['dir1'], w.genome),
        npieces = config['wgc']['npieces'],
    shell:
        """
        bed.py filter -min 5000 {params.odir}/15_intervals/11.gap.bed \
                > {params.odir}/81.qry.gap.bed
        subtractBed -nonamecheck -a {params.odir}/15_intervals/01.chrom.bed \
                -b {params.odir}/81.qry.gap.bed | bed filter -min 100 - | \
                bed makewindow -w 1000000 -s 995000 - \
                > {params.odir}/85.qry.clean.bed
        bed.py size {params.odir}/85.qry.clean.bed
        faSplit.R {params.odir}/85.qry.clean.bed {params.odir}/10_genome.fna \
                {params.odir}/15_intervals/01.chrom.sizes \
                {params.odir}/87_qrys \
                --chain {params.odir}/86.chain -n {params.npieces}
        """

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
        e = lambda w: "%s/%s/%s.%s/%s.%s.e" % (config['dirp'], config['wgc']['align']['id'], w.qry, w.tgt, w.idx, w.tchrom),
        o = lambda w: "%s/%s/%s.%s/%s.%s.o" % (config['dirp'], config['wgc']['align']['id'], w.qry, w.tgt, w.idx, w.tchrom),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'align')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'align')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'align')['mem']
    threads: config['wgc']['align']['ppn']
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
                tchrom = config['B73']['chroms'].split())
    output:
        "%s/{qry}_{tgt}/02.coord.psl" % config['wgc']['dir3']
    params:
        odir = lambda w: "%s/%s_%s" % (config['wgc']['dir3'], w.qry, w.tgt),
        chain = lambda w: "%s/%s/86.chain" %
            (config['wgc']['dir1'], w.qry),
        tbit = lambda w: config[w.tgt]['blat']['x.2bit'],
        qbit = lambda w: config[w.qry]['blat']['x.2bit'],
        tsize = lambda w: config[w.tgt]['size'],
        qsize = lambda w: config[w.qry]['size'],
        N = lambda w: "%s.%s.%s" % (config['wgc']['merge']['id'], w.qry, w.tgt),
        e = lambda w: "%s/%s/%s.%s.e" % (config['dirp'], config['wgc']['merge']['id'], w.qry, w.tgt),
        o = lambda w: "%s/%s/%s.%s.o" % (config['dirp'], config['wgc']['merge']['id'], w.qry, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'merge')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'merge')['mem']
    threads: config['wgc']['merge']['ppn']
    shell:
        """
        mkdir -p {params.odir}
        pslCat -nohead {input} > {params.odir}/01.psl
        psl.py coordQ {params.odir}/01.psl {params.qsize} \
                >{params.odir}/02.coord.psl
        pslCheck {params.odir}/02.coord.psl
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
        "%s/{qry}_{tgt}/02.coord.psl" % config['wgc']['dir3']
    output:
        "%s/{qry}_{tgt}/15.chain" % config['wgc']['dir3']
    params:
        odir = lambda w: "%s/%s_%s" % (config['wgc']['dir3'], w.qry, w.tgt),
        tbit = lambda w: config[w.tgt]['blat']['x.2bit'],
        qbit = lambda w: config[w.qry]['blat']['x.2bit'],
        tsize = lambda w: config[w.tgt]['size'],
        qsize = lambda w: config[w.qry]['size'],
        N = lambda w: "%s.%s.%s" % (config['wgc']['chain']['id'], w.qry, w.tgt),
        e = lambda w: "%s/%s/%s.%s.e" % (config['dirp'], config['wgc']['chain']['id'], w.qry, w.tgt),
        o = lambda w: "%s/%s/%s.%s.o" % (config['dirp'], config['wgc']['chain']['id'], w.qry, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'chain')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'chain')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'chain')['mem']
    threads: config['wgc']['chain']['ppn']
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

rule wgc5_post:
    input:
        "%s/{qry}_{tgt}/15.chain" % config['wgc']['dir3']
    output:
        "%s/{qry}_{tgt}/10.vnt.bed" % (config['dird']),
        "%s/{qry}_{tgt}/10.{tgt}.vcf.gz" % (config['dird']),
        "%s/{qry}_{tgt}/10.{qry}.vcf.gz" % (config['dird']),
    params:
        odir = lambda w: "%s/%s_%s" % (config['dird'], w.qry, w.tgt),
        tfas = lambda w: config[w.tgt]['ref'],
        qfas = lambda w: config[w.qry]['ref'],
        tsize = lambda w: config[w.tgt]['size'],
        qsize = lambda w: config[w.qry]['size'],
        tref = lambda w: config[w.tgt]['gatk']['xref'],
        qref = lambda w: config[w.qry]['gatk']['xref'],
        N = lambda w: "%s.%s.%s" % (config['wgc']['chain']['id'], w.qry, w.tgt),
        e = lambda w: "%s/%s/%s.%s.e" % (config['dirp'], config['wgc']['post']['id'], w.qry, w.tgt),
        o = lambda w: "%s/%s/%s.%s.o" % (config['dirp'], config['wgc']['post']['id'], w.qry, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'wgc', 'post')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'wgc', 'post')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'wgc', 'post')['mem']
    threads: config['wgc']['post']['ppn']
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
                {params.odir}/10.{wildcards.tgt}.1.vcf \
                {params.odir}/10.{wildcards.qry}.1.vcf
        gatk UpdateVCFSequenceDictionary -V \
                {params.odir}/10.{wildcards.tgt}.1.vcf -R {params.tref} \
                -O {params.odir}/10.{wildcards.tgt}.2.vcf
        gatk UpdateVCFSequenceDictionary -V \
                {params.odir}/10.{wildcards.qry}.1.vcf -R {params.qref} \
                -O {params.odir}/10.{wildcards.qry}.2.vcf
        bcftools norm -f {params.tfas} -c w -d all \
                {params.odir}/10.{wildcards.tgt}.2.vcf -Ou |\
                bcftools sort -Oz
                -o {params.odir}/10.{wildcards.tgt}.vcf.gz
        bcftools norm -f {params.qfas} -c w -d all \
                {params.odir}/10.{wildcards.qry}.2.vcf -Ou |\
                bcftools sort -Oz
                -o {params.odir}/10.{wildcards.qry}.vcf.gz
        rm {params.odir}/10.{wildcards.tgt}.[12].*
        rm {params.odir}/10.{wildcards.qry}.[12].*
        """

rule wgc6_eff_t:
    input:
        "%s/{qry}_{tgt}/10.{tgt}.vcf.gz" % config['dird'],
    output:
        "%s/{qry}_{tgt}/15.{tgt}.tsv" % config['dird'],
    params:
        odir = lambda w: "%s/%s_%s" % (config['dird'], w.qry, w.tgt),
        tfas = lambda w: config[w.tgt]['ref'],
        qfas = lambda w: config[w.qry]['ref'],
        tsize = lambda w: config[w.tgt]['size'],
        qsize = lambda w: config[w.qry]['size'],
        tref = lambda w: config[w.tgt]['gatk']['xref'],
        qref = lambda w: config[w.qry]['gatk']['xref'],
        txcfg = lambda w: config[w.tgt]['snpeff'],
        N = lambda w: "%s.%s.%s.t" % (config['wgc']['eff']['id'], w.qry, w.tgt),
        e = lambda w: "%s/%s/%s.%s.t.e" % (config['dirp'], config['wgc']['eff']['id'], w.qry, w.tgt),
        o = lambda w: "%s/%s/%s.%s.t.o" % (config['dirp'], config['wgc']['eff']['id'], w.qry, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'snpeff')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'snpeff')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'snpeff')['mem']
    threads: config['snpeff']['ppn']
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.txcfg} {wildcards.tgt} -ud 0 \
                {input[0]} | bgzip > {params.odir}/15.{wildcards.tgt}.vcf.gz
        wgc.py parseEff {params.odir}/15.{wildcards.tgt}.vcf.gz > {output[0]}
        """

rule wgc6_eff_q:
    input:
        "%s/{qry}_{tgt}/10.{qry}.vcf.gz" % config['dird'],
    output:
        "%s/{qry}_{tgt}/15.{qry}.tsv" % config['dird'],
    params:
        odir = lambda w: "%s/%s_%s" % (config['dird'], w.qry, w.tgt),
        tfas = lambda w: config[w.tgt]['ref'],
        qfas = lambda w: config[w.qry]['ref'],
        tsize = lambda w: config[w.tgt]['size'],
        qsize = lambda w: config[w.qry]['size'],
        tref = lambda w: config[w.tgt]['gatk']['xref'],
        qref = lambda w: config[w.qry]['gatk']['xref'],
        qxcfg = lambda w: config[w.qry]['snpeff'],
        N = lambda w: "%s.%s.%s.q" % (config['wgc']['eff']['id'], w.qry, w.tgt),
        e = lambda w: "%s/%s/%s.%s.q.e" % (config['dirp'], config['wgc']['eff']['id'], w.qry, w.tgt),
        o = lambda w: "%s/%s/%s.%s.q.o" % (config['dirp'], config['wgc']['eff']['id'], w.qry, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'snpeff')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'snpeff')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'snpeff')['mem']
    threads: config['snpeff']['ppn']
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.qxcfg} {wildcards.qry} -ud 0 \
                {input[0]} | bgzip > {params.odir}/15.{wildcards.qry}.vcf.gz
        wgc.py parseEff {params.odir}/15.{wildcards.qry}.vcf.gz > {output[0]}
        """

rule wgc_orthofinder:
    input:
        tgt_faa = lambda w: config[w.tgt]['lfaa'],
        qry_faa = lambda w: config[w.qry]['lfaa'],
    output:
        "%s/{qry}_{tgt}/10.tsv" % config['wgc']['dir4'],
    params:
        odir = lambda w: "%s/%s_%s" % (config['wgc']['dir4'], w.qry, w.tgt),
        odir1 = lambda w: "%s/%s_%s/01_seqs" % (config['wgc']['dir4'], w.qry, w.tgt),
        tfaa = lambda w: "%s/%s_%s/%s.faa" % (config['wgc']['dir4'], w.qry, w.tgt, w.tgt),
        qfaa = lambda w: "%s/%s_%s/%s.faa" % (config['wgc']['dir4'], w.qry, w.tgt, w.qry),
        N = lambda w: "%s.%s.%s" % (config['orthofinder']['id'], w.qry, w.tgt),
        e = lambda w: "%s/%s/%s.%s.e" % (config['dirp'], config['orthofinder']['id'], w.qry, w.tgt),
        o = lambda w: "%s/%s/%s.%s.o" % (config['dirp'], config['orthofinder']['id'], w.qry, w.tgt),
        ppn = lambda w, resources: resources.ppn,
        runtime = lambda w, resources: resources.runtime,
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config, attempt, 'orthofinder')['ppn'],
        runtime = lambda w, attempt: get_resource(config, attempt, 'orthofinder')['runtime'],
        mem = lambda w, attempt: get_resource(config, attempt, 'orthofinder')['mem']
    threads: config['orthofinder']['ppn']
    shell:
        """
        mkdir -p {params.odir}
        cp -f {input.tgt_faa} {params.tfaa}
        cp -f {input.qry_faa} {params.qfaa}
        
        source activate py27
        orthofinder -f {params.odir} -t {threads} -p {params.odir} -og
        """

