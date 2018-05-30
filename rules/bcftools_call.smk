rule bcftools_call:
    input:
        ref=config['bcftools_call']['ref'],
        samples=expand(["22.bam/{sid}_{gt}.bam"], zip, 
                sid = config['t']['sid'],
                gt = config['t']['genotype']),
        indexes=expand(["22.bam/{sid}_{gt}.bam.bai"], zip, 
                sid = config['t']['sid'],
                gt = config['t']['genotype']),
    output:
        "25.vcf"
    params:
        mpileup="",
        call=""
    log:
        "%s/bcftools_call.log" % config['dirl']
    wrapper:
        "0.23.1/bio/bcftools/call"
