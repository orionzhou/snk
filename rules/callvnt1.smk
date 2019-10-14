include: 'gatk.smk'

def cv11_inputs(w):
    yid, gt = w.yid, w.gt
    sids = config['y'][yid]['gt'][gt]
    return expand(ancient("%s/%s/{sid}.bam" % (yid, config['callvnt']['idir'])), sid = sids)

rule cv11_hc:
    input: cv11_inputs
    output:
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz" % config['callvnt']['od25']),
        temp("{yid}/%s/{gt}/{rid}.g.vcf.gz.tbi" % config['callvnt']['od25'])
    params:
        cmd = config['gatk']['cmd'],
        ref = lambda w: config['g'][config['y'][w.yid]['ref']]['gatk']['xref'],
        input_str = lambda w, input: ["-I %s" % x for x in input],
        region = lambda w: config['g'][config['y'][w.yid]['ref']]['win56'][w.rid],
        extra = gatk_extra(picard = False, jdk = True, hc = True),
        N = "{yid}.%s.{gt}.{rid}" % config['cv11_hc']['id'],
        e = "{yid}/%s/%s/{gt}/{rid}.e" % (config['dirj'], config['cv11_hc']['id']),
        o = "{yid}/%s/%s/{gt}/{rid}.o" % (config['dirj'], config['cv11_hc']['id']),
        j = lambda w: get_resource(w, config, 'cv11_hc'),
        mem = lambda w: get_resource(w, config, 'cv11_hc')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv11_hc')['ppn']
    conda: "../envs/work.yml"
    shell:
        #-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" HaplotypeCaller \
        {params.extra} -ERC GVCF \
        -R {params.ref} \
        -L {params.region} \
        {params.input_str} -O {output[0]}
        """

def cv12_inputs(w):
    yid, gt = w.yid, w.gt
    rids = natsorted(config['g'][config['y'][yid]['ref']]['win56'].keys())
    vcfs = expand("%s/%s/%s/{rid}.g.vcf.gz" %
        (yid, config['callvnt']['od25'], gt), rid = rids)
    tbis = expand("%s/%s/%s/{rid}.g.vcf.gz.tbi" %
        (yid, config['callvnt']['od25'], gt), rid = rids)
    return {'vcfs':vcfs, 'tbis':tbis}

rule cv12_merge_vcfs:
    input: unpack(cv12_inputs)
    output:
        vcf = "{yid}/%s/{gt}.g.vcf.gz" % config['callvnt']['od26'],
        tbi = "{yid}/%s/{gt}.g.vcf.gz.tbi" % config['callvnt']['od26']
    params:
        cmd = config['gatk']['cmd'],
        input_str = lambda w, input: ["-I %s" % x for x in input.vcfs],
        extra = gatk_extra(picard = True),
        N = "{yid}.%s.{gt}" % config['cv12_merge_vcfs']['id'],
        e = "{yid}/%s/%s/{gt}.e" % (config['dirj'], config['cv12_merge_vcfs']['id']),
        o = "{yid}/%s/%s/{gt}.o" % (config['dirj'], config['cv12_merge_vcfs']['id']),
        j = lambda w: get_resource(w, config, 'cv12_merge_vcfs'),
        mem = lambda w: get_resource(w, config, 'cv12_merge_vcfs')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv12_merge_vcfs')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" MergeVcfs \
        {params.extra} \
        {params.input_str} -O {output.vcf}
        """

rule cv13_rename:
    input:
        vcf = "{yid}/%s/{gt}.g.vcf.gz" % config['callvnt']['od26'],
        tbi = "{yid}/%s/{gt}.g.vcf.gz.tbi" % config['callvnt']['od26']
    output:
        vcf = protected("%s/{yid}/{gt}.g.vcf.gz" % config['callvnt']['odir']),
        tbi = protected("%s/{yid}/{gt}.g.vcf.gz.tbi" % config['callvnt']['odir']),
        stat = protected("%s/{yid}/{gt}.txt" % config['callvnt']['odir'])
    params:
        cmd = config['gatk']['cmd'],
        extra = gatk_extra(picard = True, jdk = False),
        ngt = "{yid}#{gt}",
        N = "{yid}.%s.{gt}" % config['cv13_rename']['id'],
        e = "{yid}/%s/%s/{gt}.e" % (config['dirj'], config['cv13_rename']['id']),
        o = "{yid}/%s/%s/{gt}.o" % (config['dirj'], config['cv13_rename']['id']),
        j = lambda w: get_resource(w, config, 'cv13_rename'),
        mem = lambda w: get_resource(w, config, 'cv13_rename')['mem'],
    resources: attempt = lambda w, attempt: attempt
    threads: lambda w: get_resource(w, config, 'cv13_rename')['ppn']
    conda: "../envs/work.yml"
    shell:
        """
        {params.cmd} --java-options "-Xmx{params.mem}G" RenameSampleInVcf \
            {params.extra} \
            --NEW_SAMPLE_NAME {params.ngt} --OLD_SAMPLE_NAME {wildcards.gt} \
            -I {input.vcf} -O {output.vcf} --CREATE_INDEX
        bcftools stats -s - {output.vcf} > {output.stat}
        """

rule cv14_vcfstats:
    input:
        lambda w: expand("%s/%s/{gt}.txt" % (config['callvnt']['odir'], w.yid), gt = config['y'][w.yid]['Genotypes'])
    output:
        snp = "%s/{yid}/%s" % (config['oid'], config['callvnt']['out']),
        idl = "%s/{yid}/%s" % (config['oid'], config['callvnt']['out_idl'])
    shell: """
        grep '^PSC' {input} > {output.snp}
        grep '^PSI' {input} > {output.idl}"""



