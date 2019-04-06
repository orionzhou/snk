rule bcftools_call:
    input:
        lambda w: expand("%s/%s/{sid}.bam" % (w.yid, config['cleanbam']['od23']),
            sid = config['y'][w.yid]['SampleID'])
    output:
        protected("{yid}/%s/{rid}.vcf.gz" % config['cleanbam']['od24a']),
        protected("{yid}/%s/{rid}.vcf.gz.tbi" % config['cleanbam']['od24a'])
    params:
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['fasta']['ref'],
        region = lambda w: config['g'][config['y'][w.yid]['reference']]['regions'][w.rid],
        N = "{yid}.%s.{rid}" % config['bcftools']['call']['id'],
        e = "{yid}/%s/%s/{rid}.e" % (config['dirj'], config['bcftools']['call']['id']),
        o = "{yid}/%s/%s/{rid}.o" % (config['dirj'], config['bcftools']['call']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'call')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'call')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'call')['mem']
    threads: config["bcftools"]['call']["ppn"]
    shell:
        """
        bcftools mpileup -f {params.ref} -r {params.region} \
            {input} -Ou | \
            bcftools call -c -Oz -o {output[0]} -
        bcftools index -t {output[0]}
        """

rule bcftools_concat:
    input:
        lambda w: expand("%s/%s/{rid}.vcf.gz" % (w.yid, config['cleanbam']['od24a']),
            rid = config['g'][config['y'][w.yid]['reference']]['regions'].keys())
    output:
        vcf = protected("{yid}/%s" % config['cleanbam']['of24b']),
        tbi = protected("{yid}/%s.tbi" % config['cleanbam']['of24b']),
        stat = protected("{yid}/%s" % config['cleanbam']['of24b'].replace(".vcf.gz", ".txt"))
    params:
        ref = lambda w: config['g'][config['y'][w.yid]['reference']]['fasta']['ref'],
        N = "{yid}.%s" % config['bcftools']['concat']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['bcftools']['concat']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['bcftools']['concat']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'concat')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'concat')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'concat')['mem']
    threads:
        config["bcftools"]['concat']["ppn"]
    shell:
        """
        bcftools concat {input} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        bcftools stats -s - {output.vcf} > {output.stat}
        """

rule bcftools_filter:
    input:
        vcf = "{yid}/%s" % config['cleanbam']['of24b'],
        tbi = "{yid}/%s.tbi" % config['cleanbam']['of24b']
    output:
        vcf = protected("{yid}/%s" % config['cleanbam']['of24c']),
        tbi = protected("{yid}/%s.tbi" % config['cleanbam']['of24c']),
        stat = protected("{yid}/%s" % config['cleanbam']['of24c'].replace(".vcf.gz", ".txt"))
    params:
        minqual = 900,
        N = "{yid}.%s" % config['bcftools']['filter']['id'],
        e = "{yid}/%s/%s.e" % (config['dirj'], config['bcftools']['filter']['id']),
        o = "{yid}/%s/%s.o" % (config['dirj'], config['bcftools']['filter']['id']),
    resources:
        ppn = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'filter')['ppn'],
        runtime = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'filter')['runtime'],
        mem = lambda w, attempt:  get_resource(config, attempt, 'bcftools', 'filter')['mem']
    threads:
        config["bcftools"]['filter']["ppn"]
    shell:
        """
        bcftools view -i 'QUAL>{params.minqual}' {input.vcf} -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        bcftools stats -s - {output.vcf} > {output.stat}
        """

rule bcftools_gtcheck:
    input:
        "%s/{sid}.bcf" % config['cleanbam']['od23']
    output:
        "%s/{sid}.txt" % config['cleanbam']['od23']
    params:
        vcf = '',
        extra = '',
    shell:
        """
        bcftools gtcheck {params.extra} -g {params.vcf} {input} > {output}
        """

