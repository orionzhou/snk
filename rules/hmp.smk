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
        ppn = lambda w, attempt: get_resource(w, config,attempt,'hmp','crossmap')['ppn'],
        runtime = lambda w, attempt: get_resource(w, config,attempt,'hmp','crossmap')['runtime'],
        mem = lambda w, attempt: get_resource(w, config,attempt,'hmp','crossmap')['mem']
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
        ppn = lambda w, attempt: get_resource(w, config,attempt,'hmp','vcfnorm')['ppn'],
        runtime = lambda w, attempt: get_resource(w, config,attempt,'hmp','vcfnorm')['runtime'],
        mem = lambda w, attempt: get_resource(w, config,attempt,'hmp','vcfnorm')['mem']
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
        ppn = lambda w, attempt: get_resource(w, config,attempt,'hmp','liftover')['ppn'],
        runtime = lambda w, attempt: get_resource(w, config,attempt,'hmp','liftover')['runtime'],
        mem = lambda w, attempt: get_resource(w, config,attempt,'hmp','liftover')['mem']
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
        ppn = lambda w, attempt: get_resource(w, config,attempt,'hmp','merge')['ppn'],
        runtime = lambda w, attempt: get_resource(w, config,attempt,'hmp','merge')['runtime'],
        mem = lambda w, attempt: get_resource(w, config,attempt,'hmp','merge')['mem']
    threads: config['hmp']['merge']['ppn']
    conda: "envs/samtools.yml"
    shell:
        """
        bcftools concat --naive -Oz -o {params.o1} {input}
        bcftools sort -m {params.mem}G -T {params.tmp} -Oz -o {output} {params.o1}
        bcftools index -t {output}
        """


