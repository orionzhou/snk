rule wgc_prep_tgt:
    input: fna = lambda w: config['g'][w.genome]['fasta']['ref'],
    output: "%s/{genome}/11_chroms/{tchrom}.fna" % config['wgc']['od10']
    params:
        odir = "%s/{genome}" % config['wgc']['od10'],
        cdir = "%s/{genome}/11_chroms" % config['wgc']['od10'],
        N = "%s.{genome}" % config['wgc_prep_tgt']['id'],
        e = "%s/%s/{genome}.e" % (config['dirj'], config['wgc_prep_tgt']['id']),
        o = "%s/%s/{genome}.o" % (config['dirj'], config['wgc_prep_tgt']['id']),
        j = lambda w: get_resource(w, config, 'wgc_prep_tgt'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_prep_tgt']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        mkdir -p {params.cdir}
        fasta.py extract {input.fna} {wildcards.tchrom} > {output}
        """

rule wgc_prep_qry:
    input:
        fna = lambda w: config['g'][w.genome]['fasta']['ref'],
        chrom_size = lambda w: config['g'][w.genome]['fasta']['chrom_size'],
        chrom_bed = lambda w: config['g'][w.genome]['fasta']['chrom_bed'],
        gap = lambda w: config['g'][w.genome]['fasta']['gap'],
    output:
        fnas = expand("%s/{{genome}}/87_qrys/part.{idx}.fna" % \
                config['wgc']['od10'],
                idx = range(1,config['wgc']['npieces']+1)),
        chain = "%s/{genome}/86.chain" % config['wgc']['od10'],
    params:
        odir = "%s/{genome}" % config['wgc']['od10'],
        cdir = "%s/{genome}/11_chroms" % config['wgc']['od10'],
        npieces = config['wgc']['npieces'],
        N = "%s.{genome}" % config['wgc_prep_qry']['id'],
        e = "%s/%s/{genome}.e" % (config['dirj'], config['wgc_prep_qry']['id']),
        o = "%s/%s/{genome}.o" % (config['dirj'], config['wgc_prep_qry']['id']),
        j = lambda w: get_resource(w, config, 'wgc_prep_qry'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_prep_qry']['ppn']
    conda: "../envs/work.yml"
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

rule wgc_align:
    input:
        tgt = "%s/{tgt}/11_chroms/{tchrom}.fna" % config['wgc']['od10'],
        qry = "%s/{qry}/87_qrys/part.{idx}.fna" % config['wgc']['od10'],
    output:
        "%s/{qry}_{tgt}/q{idx}.{tchrom}.psl" % config['wgc']['od20']
    params:
        sam = "%s/{qry}_{tgt}/q{idx}.{tchrom}.sam" % config['wgc']['od20'],
        lav = "%s/{qry}_{tgt}/q{idx}.{tchrom}.lav" % config['wgc']['od20'],
        N = lambda w: "%s.%s%s.%s.%s" % (config['wgc_align']['id'], w.qry[:1], w.tgt[:1], w.idx, w.tchrom),
        e = lambda w: "%s/%s/%s.%s/%s.%s.e" % (config['dirj'], config['wgc_align']['id'], w.qry, w.tgt, w.idx, w.tchrom),
        o = lambda w: "%s/%s/%s.%s/%s.%s.o" % (config['dirj'], config['wgc_align']['id'], w.qry, w.tgt, w.idx, w.tchrom),
        j = lambda w: get_resource(w, config, 'wgc_align')
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_align']['ppn']
    conda: "../envs/work.yml"
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

def wgc_merge_inputs(w):
    odir = "%s/%s_%s" % (config['wgc']['od20'], w.qry, w.tgt)
    fi = config['g'][w.tgt]['fasta']['chrom_size']
    df = pd.read_csv(fi, sep='\t', names=['chrom','size'])
    tchroms = df['chrom'].values.tolist()
    return expand("%s/q{idx}.{tchrom}.psl" % odir,
                  idx = range(1, config['wgc']['npieces']+1),
                  tchrom = tchroms)

rule wgc_merge:
    input: wgc_merge_inputs
    output: "%s/{qry}_{tgt}/03.coord.pass.psl" % config['wgc']['od23']
    params:
        odir = "%s/{qry}_{tgt}" % config['wgc']['od23'],
        chain = lambda w: "%s/%s/86.chain" %
            (config['wgc']['od10'], w.qry),
        tbit = lambda w: config['g'][w.tgt]['blat']['x.2bit'],
        qbit = lambda w: config['g'][w.qry]['blat']['x.2bit'],
        tsize = lambda w: config['g'][w.tgt]['fasta']['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['fasta']['chrom_size'],
        N = "%s.{qry}.{tgt}" % config['wgc_merge']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc_merge']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc_merge']['id']),
        j = lambda w: get_resource(w, config, 'wgc_prep_tgt'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_merge']['ppn']
    conda: "../envs/work.yml"
    shell:
#        pslSwap {params.odir}/01.psl {params.odir}/02.swap.psl
        """
        pslCat -nohead {input} > {params.odir}/01.psl
        psl.py coordQ {params.odir}/01.psl {params.qsize} \
                >{params.odir}/02.coord.psl
        mkdir -p {params.odir}
        pslCheck {params.odir}/02.coord.psl \
            -pass={params.odir}/02.coord.pass.psl \
            -fail={params.odir}/02.coord.fail.psl || echo non_success

        psl.py coordT {params.odir}/02.coord.psl {params.tsize} \
                >{params.odir}/03.coord.psl
        pslCheck {params.odir}/03.coord.psl \
                -pass={params.odir}/03.coord.pass.psl \
                -fail={params.odir}/03.coord.fail.psl || echo non_success
        """

rule wgc_chain:
    input: "%s/{qry}_{tgt}/03.coord.pass.psl" % config['wgc']['od23']
    output: "%s/{qry}_{tgt}/15.chain" % config['wgc']['od23']
    params:
        odir = "%s/{qry}_{tgt}" % config['wgc']['od23'],
        tbit = lambda w: config['g'][w.tgt]['blat']['x.2bit'],
        qbit = lambda w: config['g'][w.qry]['blat']['x.2bit'],
        tsize = lambda w: config['g'][w.tgt]['fasta']['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['fasta']['chrom_size'],
        N = "%s.{qry}.{tgt}" % config['wgc_chain']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc_chain']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc_chain']['id']),
        j = lambda w: get_resource(w, config, 'wgc_chain'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_chain']['ppn']
    conda: "../envs/work.yml"
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

rule wgc_post: # python+gatk+R
    input: "%s/{qry}_{tgt}/15.chain" % config['wgc']['od23']
    output:
        "%s/{qry}_{tgt}/10.vnt.bed" % config['dirr'],
        "%s/{qry}_{tgt}/10.{tgt}.vcf.gz" % config['dirr'],
        "%s/{qry}_{tgt}/10.{qry}.vcf.gz" % config['dirr'],
    params:
        odir = "%s/{qry}_{tgt}" % config['dirr'],
        tpre = "%s/{qry}_{tgt}/10.{tgt}" % config['dirr'],
        qpre = "%s/{qry}_{tgt}/10.{qry}" % config['dirr'],
        tfas = lambda w: config['g'][w.tgt]['fasta']['ref'],
        qfas = lambda w: config['g'][w.qry]['fasta']['ref'],
        tsize = lambda w: config['g'][w.tgt]['fasta']['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['fasta']['chrom_size'],
        tref = lambda w: config['g'][w.tgt]['gatk']['xref'],
        qref = lambda w: config['g'][w.qry]['gatk']['xref'],
        N = "%s.{qry}.{tgt}" % config['wgc_post']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc_post']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc_post']['id']),
        j = lambda w: get_resource(w, config, 'wgc_post'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_post']['ppn']
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

rule wgc_eff_t: # python+snpeff
    input: "%s/{qry}_{tgt}/10.{tgt}.vcf.gz" % config['dirr'],
    output: "%s/{qry}_{tgt}/15.{tgt}.tsv" % config['dirr'],
    params:
        odir = "%s/{qry}_{tgt}" % config['dirr'],
        tfas = lambda w: config['g'][w.tgt]['fasta']['ref'],
        qfas = lambda w: config['g'][w.qry]['fasta']['ref'],
        tsize = lambda w: config['g'][w.tgt]['fasta']['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['fasta']['chrom_size'],
        tref = lambda w: config['g'][w.tgt]['gatk']['xref'],
        qref = lambda w: config['g'][w.qry]['gatk']['xref'],
        txcfg = lambda w: config['g'][w.tgt]['snpeff']['xcfg'],
        N = "%s.{qry}.{tgt}.t" % config['wgc_eff']['id'],
        e = "%s/%s/{qry}.{tgt}.t.e" % (config['dirj'], config['wgc_eff']['id']),
        o = "%s/%s/{qry}.{tgt}.t.o" % (config['dirj'], config['wgc_eff']['id']),
        j = lambda w: get_resource(w, config, 'wgc_eff'),
        mem = lambda w: get_resource(w, config, 'wgc_eff')['mem'],
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_eff']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.txcfg} {wildcards.tgt} -ud 0 \
                {input[0]} | bgzip > {params.odir}/15.{wildcards.tgt}.vcf.gz
        wgc.py parseEff {params.odir}/15.{wildcards.tgt}.vcf.gz > {output[0]}
        """

rule wgc_eff_q:
    input: "%s/{qry}_{tgt}/10.{qry}.vcf.gz" % config['dirr'],
    output: "%s/{qry}_{tgt}/15.{qry}.tsv" % config['dirr'],
    params:
        odir = "%s/{qry}_{tgt}" % config['dirr'],
        tfas = lambda w: config['g'][w.tgt]['fasta']['ref'],
        qfas = lambda w: config['g'][w.qry]['fasta']['ref'],
        tsize = lambda w: config['g'][w.tgt]['fasta']['chrom_size'],
        qsize = lambda w: config['g'][w.qry]['fasta']['chrom_size'],
        tref = lambda w: config['g'][w.tgt]['gatk']['xref'],
        qref = lambda w: config['g'][w.qry]['gatk']['xref'],
        qxcfg = lambda w: config['g'][w.qry]['snpeff']['xcfg'],
        N = "%s.{qry}.{tgt}.q" % config['wgc_eff']['id'],
        e = "%s/%s/{qry}.{tgt}.q.e" % (config['dirj'], config['wgc_eff']['id']),
        o = "%s/%s/{qry}.{tgt}.q.o" % (config['dirj'], config['wgc_eff']['id']),
        j = lambda w: get_resource(w, config, 'wgc_eff'),
        mem = lambda w: get_resource(w, config, 'wgc_eff')['mem'],
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_eff']['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        snpEff -Xmx{params.mem}G -c {params.qxcfg} {wildcards.qry} -ud 0 \
                {input[0]} | bgzip > {params.odir}/15.{wildcards.qry}.vcf.gz
        wgc.py parseEff {params.odir}/15.{wildcards.qry}.vcf.gz > {output[0]}
        """

rule wgc_blastp:
    input:
        tgt_faa = lambda w: config['g'][w.tgt]['annotation']['lfaa'],
        qry_faa = lambda w: config['g'][w.qry]['annotation']['lfaa'],
        tgt_db = lambda w: config['g'][w.tgt]['blastp']['xout'],
    output: "%s/{qry}_{tgt}.tsv" % config['wgc']['od42'],
    params:
        tgt_db = lambda w: config['g'][w.tgt]['blastp']['xpre'],
        N = "%s.{qry}.{tgt}" % config['wgc_blastp']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc_blastp']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc_blastp']['id']),
        j = lambda w: get_resource(w, config, 'wgc_blastp'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_blastp']['ppn']
    conda: "../envs/blast.yml"
    shell:
        """
        blastp -db {params.tgt_db} -query {input.qry_faa} -num_threads {threads} -outfmt 6 -out {output}
        """

rule wgc_synteny:
    input:
        tgt_faa = lambda w: config['g'][w.tgt]['annotation']['lfaa'],
        qry_faa = lambda w: config['g'][w.qry]['annotation']['lfaa'],
        tgt_bed = lambda w: config['g'][w.tgt]['annotation']['lbed'],
        qry_bed = lambda w: config['g'][w.qry]['annotation']['lbed'],
        last = "%s/{qry}_{tgt}.tsv" % config['wgc']['od42'],
    output: "%s/{qry}_{tgt}/%s/%s" % (config['dirr'], config['wgc']['rd20'], config['wgc']['out']),
    params:
        in_last = lambda w, input: op.abspath(input.last),
        out = lambda w, output: op.basename(output[0]),
        odir = "%s/{qry}_{tgt}/%s" % (config['dirr'], config['wgc']['rd20']),
        cscore = 0.7,
        quota = '1:1',
        dist = 20,
        Nm = 10,
        N = "%s.{qry}.{tgt}" % config['wgc_synteny']['id'],
        e = "%s/%s/{qry}.{tgt}.e" % (config['dirj'], config['wgc_synteny']['id']),
        o = "%s/%s/{qry}.{tgt}.o" % (config['dirj'], config['wgc_synteny']['id']),
        j = lambda w: get_resource(w, config, 'wgc_synteny'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['wgc_synteny']['ppn']
    conda: "../envs/work.yml"
    script: "../scripts/synteny.py"

rule orthofinder:
    input: [config['g'][genome]['annotation']['lfaa'] for genome in config['ortho_genomes']]
    output: "%s/%s" % (config['wgc']['od50'], config['wgc']['of50'])
    params:
        odir = config['wgc']['od50'],
        N = "%s" % config['orthofinder']['id'],
        e = "%s/%s.e" % (config['dirj'], config['orthofinder']['id']),
        o = "%s/%s.o" % (config['dirj'], config['orthofinder']['id']),
        j = lambda w: get_resource(w, config, 'orthofinder'),
    resources:
        attempt = lambda w, attempt: attempt,
    threads: config['orthofinder']['ppn']
    conda: "../envs/orthofinder.yml"
    script: "../scripts/run_orthofinder.py"


