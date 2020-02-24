rule brb1_trim:
    input:
        fq1 = "%s/{sid}_1.fq.gz" % config['brbseq']['idir'],
        fq2 = "%s/{sid}_2.fq.gz" % config['brbseq']['idir']
    output:
        fq1 = "%s/{sid}_1.fq.gz" % config['brbseq']['od02'],
        fq2 = "%s/{sid}_2.fq.gz" % config['brbseq']['od02']
    params:
        brb = 'java -jar /home/springer/zhoux379/source/git/BRB-seqTools/releases/BRBseqTools.1.4.jar',
        len1 = 21,
        umi = 15,
        od02 = config['brbseq']['od02'],
        trim2 = "%s/{sid}_2.trimmed.fastq.gz" % config['brbseq']['od02'],
        N = "%s.{sid}" % config['brb1_trim']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['brb1_trim']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['brb1_trim']['id']),
        j = lambda w: get_resource(w, config, 'brb1_trim'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'brb1_trim')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        fastp --thread {threads} --disable_quality_filtering \
            --disable_trim_poly_g --max_len1 {params.len1} \
            -i {input.fq1} -o {output.fq1}

        {params.brb} Trim -f {input.fq2} -o {params.od02}
        mv {params.trim2} {output.fq2}

        fastqc --threads {threads} --noextract --format fastq -o {params.od02} {input.fq2}
        """

rule brb1b_umi:
    input:
        fq1 = "%s/{sid}_1.fq.gz" % config['brbseq']['od02'],
        fq2 = "%s/{sid}_2.fq.gz" % config['brbseq']['idir']
    output:
        stat = "%s/{sid}/stats.txt" % config['brbseq']['od03'],
        cnt = "%s/{sid}/A1.tsv" % config['brbseq']['od04']
    params:
        brb = 'java -jar /home/springer/zhoux379/source/git/BRB-seqTools/releases/BRBseqTools.1.4.jar',
        len1 = 21,
        umi = 15,
        barcode = "$rn/misc/barcode1.tsv",
        od1 = "%s/{sid}" % config['brbseq']['od03'],
        od2 = "%s/{sid}" % config['brbseq']['od04'],
        N = "%s.{sid}" % config['brb1_trim']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['brb1_trim']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['brb1_trim']['id']),
        j = lambda w: get_resource(w, config, 'brb1_trim'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'brb1_trim')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        mkdir -p {params.od1} {params.od2}
        {params.brb} Demultiplex -p BU -UMI {params.umi} \
            -r1 {input.fq1} -r2 {input.fq2} \
            -c {params.barcode} -o {params.od1}
        parallel -j {threads} "fastq.py UMIcount {{}} {params.od2}/{{/..}}.tsv" ::: {params.od1}/*.fastq.gz
        """

def db_index(w, db):
    ref = 'Zmays_B73'

    if db == 'star':
        return config['g'][ref]["star"]['xpre']
    elif db == 'hisat2':
        return config['g'][ref]["hisat2"]['xpre']
    elif db == 'bwa':
        return config['g'][ref]["bwa"]['xpre']
    elif db == 'bismark':
        return config['g'][ref]["bismark"]['xpre']
    elif db == 'gtf':
        return config['g'][ref]["annotation"]['gtf']
    elif db == 'gatk':
        return config['g'][ref]["gatk"]['xref']
    else:
        print('unknown db: %d' % db)
        sys.exit(1)

def hisat2_extra(w):
    extras = "--new-summary --no-softclip".split()
    extras = "--new-summary".split()
    extras.append("--rg-id %s --rg SM:%s" % (w.sid, w.sid))
    return " ".join(extras)

rule brb2_align:
    input:
        fq2 = "%s/{sid}_2.fq.gz" % config['brbseq']['od02']
    output:
        temp("%s/{sid}.bam" % config['brbseq']['od11']),
        "%s/{sid}.txt" % config['brbseq']['od11']
    params:
        index = lambda w: db_index(w, 'hisat2'),
        extra = hisat2_extra,
        N = "%s.{sid}" % config['brb2_align']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['brb2_align']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['brb2_align']['id']),
        j = lambda w: get_resource(w, config, 'brb2_align'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'brb2_align')['ppn']
    conda: "../envs/hisat2.yml"
    shell:
        """
        hisat2 {params.extra} --threads {threads} \
            -x {params.index} -U {input.fq2} \
            --summary-file {output[1]} \
            | samtools view -Sbh -o {output[0]} -
        """

rule brb2b_sort:
    input: "%s/{sid}.bam" % config['brbseq']['od11']
    output:
        "%s/{sid}.bam" % config['brbseq']['od12'],
        "%s/{sid}.bam.bai" % config['brbseq']['od12']
    params:
        pre = "%s/{sid}" % config['brbseq']['od12'],
        tmp_dir = config['tmpdir'],
        N = "%s.{sid}" % config['brb2b_sort']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['brb2b_sort']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['brb2b_sort']['id']),
        j = lambda w: get_resource(w, config, 'brb2b_sort')
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'brb2b_sort')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        sambamba sort --tmpdir={params.tmp_dir} -t {threads} -o {output[0]} {input}
        samtools index {output[0]}
        """

rule brb2c_stat:
    input:
        fq1 = "%s/{sid}_1.fq.gz" % config['brbseq']['od02'],
        bam = "%s/{sid}.bam" % config['brbseq']['od12']
    output:
        stat = "%s/{sid}.tsv" % config['brbseq']['od13'],
        unmap = "%s/{sid}.unmap.bam" % config['brbseq']['od13'],
        abam = "%s/{sid}.bam" % config['brbseq']['od13'],
    params:
        brb = 'java -jar /home/springer/zhoux379/source/git/BRB-seqTools/releases/BRBseqTools.1.4.jar',
        abam = '{sid}.bam.annotated.bam',
        barcode = "$rn/misc/barcode1.tsv",
        gtf = lambda w: db_index(w, 'gtf'),
        tmp = "%s/brb2c_stat_{sid}" % config['tmpdir'],
        len1 = 21,
        umi = 15,
        N = "%s.{sid}" % config['brb2c_stat']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['brb2c_stat']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['brb2c_stat']['id']),
        j = lambda w: get_resource(w, config, 'brb2c_stat'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'brb2c_stat')['ppn']
    conda: "../envs/work.yml"
    shell: """
        bam.py stat {input.bam} > {output.stat}
        samtools view -b -f 4 {input.bam} > {output.unmap}

        mkdir -p {params.tmp}
        {params.brb} AnnotateBAM -p BU -UMI {params.umi} \
            -f {input.fq1} -b {input.bam} \
            -c {params.barcode} -gtf {params.gtf} -t {params.tmp} -o .
        mv {params.abam} {output.abam}

        #picard MarkDuplicates \
            #BARCODE_TAG=BC READ_ONE_BARCODE_TAG=BX \
            #I=hisat2.sorted.bam.annotated.bam O=hisat2.sorted.md.bam \
            #M=hisat2.sorted.md.txt
        """

rule brb3_dge:
    input:
        fq1 = "%s/{sid}_1.fq.gz" % config['brbseq']['od02'],
        bam = "%s/{sid}.bam" % config['brbseq']['od12']
    output:
        dge = "%s/{sid}/output.dge.umis.txt" % config['brbseq']['od21']
    params:
        brb = 'java -jar /home/springer/zhoux379/source/git/BRB-seqTools/releases/BRBseqTools.1.4.jar',
        dge = "%s/{sid}" % config['brbseq']['od21'],
        barcode = "$rn/misc/barcode1.tsv",
        gtf = lambda w: db_index(w, 'gtf'),
        tmp = "%s/{sid}/tmp" % config['brbseq']['od21'],
        len1 = 21,
        umi = 15,
        N = "%s.{sid}" % config['brb3_dge']['id'],
        e = "%s/%s/{sid}.e" % (config['dirj'], config['brb3_dge']['id']),
        o = "%s/%s/{sid}.o" % (config['dirj'], config['brb3_dge']['id']),
        j = lambda w: get_resource(w, config, 'brb3_dge'),
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'brb3_dge')['ppn'] - 1
    conda: "../envs/work.yml"
    shell: """
        rm -rf {params.tmp}
        mkdir -p {params.tmp}
        samtools view -bh -q 20 -@ {threads} -o {params.tmp}/in.bam {input.bam}
        {params.brb} CreateDGEMatrix -p BU -UMI {params.umi} \
            -f {input.fq1} -b {params.tmp}/in.bam \
            -c {params.barcode} -gtf {params.gtf} -t {params.tmp} -o {params.dge}
        rm -rf {params.tmp}
        """

