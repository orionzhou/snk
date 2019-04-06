rule crossmap:
    input: "%s/merged_flt_c{cid}.imputed.vcf.gz" % config['hmp']['od01']
    output: temp("%s/{cid}.vcf" % config['hmp']['od02'])
    params:
        refo = config['hmp']['refo'],
        refn = config['hmp']['refn'],
        chain3 = config['hmp']['chain3'],
        chain4 = config['hmp']['chain4'],
        o1 = "%s/{cid}.1.vcf" % config['hmp']['od01'],
        N = "%s.{cid}" % config['hmp']['crossmap']['id'],
        e = "%s/%s/{cid}.e" % (config['dirj'], config['hmp']['crossmap']['id']),
        o = "%s/%s/{cid}.o" % (config['dirj'], config['hmp']['crossmap']['id']),
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','crossmap')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','crossmap')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','crossmap')['mem']
    threads: config['hmp']['crossmap']['ppn']
    conda: "envs/crossmap.yml"
    shell:
        "CrossMap.py vcf {params.chain3} {input} {params.refo} {output}"

rule vcfnorm:
    input: "%s/{cid}.vcf" % config['hmp']['od02']
    output:
        vcf = temp("%s/{cid}.vcf.gz" % config['hmp']['od03']),
        tbi = temp("%s/{cid}.vcf.gz.tbi" % config['hmp']['od03']),
    params:
        refo = config['hmp']['refo'],
        refn = config['hmp']['refn'],
        chain3 = config['hmp']['chain3'],
        chain4 = config['hmp']['chain4'],
        o1 = "%s/{cid}.1.vcf.gz" % config['hmp']['od03'],
        o2 = "%s/{cid}.2.vcf.gz" % config['hmp']['od03'],
        o3 = "%s/{cid}.3.vcf.gz" % config['hmp']['od03'],
        N = "%s.{cid}" % config['hmp']['vcfnorm']['id'],
        e = "%s/%s/{cid}.e" % (config['dirj'], config['hmp']['vcfnorm']['id']),
        o = "%s/%s/{cid}.o" % (config['dirj'], config['hmp']['vcfnorm']['id']),
        mem = lambda w, resources: resources.mem - 20
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','vcfnorm')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','vcfnorm')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','vcfnorm')['mem']
    threads: config['hmp']['vcfnorm']['ppn']
    conda: "envs/samtools.yml"
    shell:
#        bcftools sort -m {params.mem}G -Oz -o {params.o1} {input}
#        bcftools norm -c x -m -any -f {params.refo} {params.o1} -Oz -o {params.o2}
#        bcftools filter -e "REF==ALT[0]" {params.o2} -Oz -o {params.o3}
        """
        bcftools norm -m +any -f {params.refo} {params.o3} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """

rule liftover:
    input:
        "%s/{cid}.vcf.gz" % config['hmp']['od03'],
        "%s/{cid}.vcf.gz.tbi" % config['hmp']['od03']
    output:
        vcf = temp("%s/{cid}.vcf.gz" % config['hmp']['od04']),
        tbi = temp("%s/{cid}.vcf.gz.tbi" % config['hmp']['od04']),
    params:
        refo = config['hmp']['refo'],
        refn = config['hmp']['refn'],
        chain3 = config['hmp']['chain3'],
        chain4 = config['hmp']['chain4'],
        rej = "%s/{cid}.rej.vcf" % config['hmp']['od04'],
        N = "%s.{cid}" % config['hmp']['liftover']['id'],
        e = "%s/%s/{cid}.e" % (config['dirj'], config['hmp']['liftover']['id']),
        o = "%s/%s/{cid}.o" % (config['dirj'], config['hmp']['liftover']['id']),
        mem = lambda w, resources: resources.mem - 5
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','liftover')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','liftover')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','liftover')['mem']
    threads: config['hmp']['liftover']['ppn']
    conda: "envs/gatk.yml"
    shell:
        """
        gatk --java-options "-Xmx{params.mem}G" LiftoverVcf \
        --MAX_RECORDS_IN_RAM 5000 \
        --TMP_DIR /scratch.global/zhoux379/temp/ \
        --USE_JDK_DEFLATER --USE_JDK_INFLATER \
        -I {input[0]} -O {output.vcf} --REJECT {params.rej} \
        -C {params.chain4} -R {params.refn}
        """

rule merge:
    input: expand("%s/{cid}.vcf.gz" % config['hmp']['od04'], cid=range(1,11))
    output: config['hmp']['of10']
    params:
        tmp = config['tmpdir'],
        o1 = "10.unsorted.vcf.gz",
        N = config['hmp']['merge']['id'],
        e = "%s/%s.e" % (config['dirj'], config['hmp']['merge']['id']),
        o = "%s/%s.o" % (config['dirj'], config['hmp']['merge']['id']),
        mem = lambda w, resources: resources.mem - 15
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','merge')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','merge')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','merge')['mem']
    threads: config['hmp']['merge']['ppn']
    conda: "envs/samtools.yml"
    shell:
        """
        bcftools concat --naive -Oz -o {params.o1} {input}
        bcftools sort -m {params.mem}G -T {params.tmp} -Oz -o {output} {params.o1}
        bcftools index -t {output}
        """

rule filter:
    input: config['hmp']['of10']
    output:
        snp = config['hmp']['of21'],
        idl = config['hmp']['of22'],
    params:
        tmp = config['tmpdir'],
        snp_stat = lambda w, output: output.snp.replace(".vcf.gz", ".stats.txt"),
        idl_stat = lambda w, output: output.idl.replace(".vcf.gz", ".stats.txt"),
        snp_sites = lambda w, output: output.snp.replace(".vcf.gz", ".sites.vcf.gz"),
        idl_sites = lambda w, output: output.idl.replace(".vcf.gz", ".sites.vcf.gz"),
        snp_bed = lambda w, output: output.snp.replace(".vcf.gz", ".bed"),
        idl_bed = lambda w, output: output.idl.replace(".vcf.gz", ".bed"),
        N = config['hmp']['filter']['id'],
        e = "%s/%s.e" % (config['dirj'], config['hmp']['filter']['id']),
        o = "%s/%s.o" % (config['dirj'], config['hmp']['filter']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','filter')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','filter')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','filter')['mem']
    threads: config['hmp']['filter']['ppn']
    conda: "envs/samtools.yml"
    shell:
#        bcftools view -m2 -M2 -v indels -i 'MAF>0.05' {input} -Oz -o {output.idl}
        """
        bcftools view -m2 -M2 -v snps -i 'MAF>0.05' {input} -Oz -o {output.snp}
        bcftools view -m2 -c 4 -v indels {input} -Oz -o {output.idl}
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
    input: config['hmp']['of21']
    output: config['hmp']['of31']
    params:
        fraction = 0.001,
        refn = config['hmp']['refn'],
        tmp = config['tmpdir'],
        N = config['hmp']['sample']['id'],
        e = "%s/%s.e" % (config['dirj'], config['hmp']['sample']['id']),
        o = "%s/%s.o" % (config['dirj'], config['hmp']['sample']['id']),
        mem = lambda w, resources: resources.mem
    resources:
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','sample')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','sample')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','sample')['mem']
    threads: config['hmp']['sample']['ppn']
    conda: "envs/gatk.yml"
    shell:
        """
        gatk --java-options "-Xmx{params.mem}G" SelectVariants \
        --tmp-dir /scratch.global/zhoux379/temp/ \
        --use-jdk-deflater --use-jdk-inflater \
        -R {params.refn} \
        -V {input} -O {output} -fraction {params.fraction}
        """

rule phylo:
    input: config['hmp']['of31']
    output: config['hmp']['of35']
    params:
        phy = lambda w, output: output.replace(".iqtree", ".phy"),
        vphy = lambda w, output: output.replace(".iqtree", ".varsites.phy"),
        tmp = config['tmpdir'],
        N = config['hmp']['phylo']['id'],
        e = "%s/%s.e" % (config['dirj'], config['hmp']['phylo']['id']),
        o = "%s/%s.o" % (config['dirj'], config['hmp']['phylo']['id']),
        ppn = lambda w, resources: resources.ppn,
        mem = lambda w, resources: resources.mem
    resources:
        q = lambda w, attempt: get_resource(config,attempt,'hmp','phylo')['q'],
        ppn = lambda w, attempt: get_resource(config,attempt,'hmp','phylo')['ppn'],
        runtime = lambda w, attempt: get_resource(config,attempt,'hmp','phylo')['runtime'],
        mem = lambda w, attempt: get_resource(config,attempt,'hmp','phylo')['mem']
    threads: config['hmp']['phylo']['ppn']
    conda: "envs/job.yml"
    shell:
        """
        vcf2phylip.py -i {input} -o {params.phy}
        iqtree -s {params.vphy} -m MFP+ASC -bb 1000 -nt {params.ppn} -pre {params.prx}
        """


