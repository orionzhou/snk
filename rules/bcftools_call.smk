rule bcftools_call:
    input:
        ref=config['bcftools']['ref'],
        samples=expand(["%s/{sid}_{gt}.bam" % config['bcftools']['call']['idir']], zip, 
                sid = config['t']['SampleID'],
                gt = config['t']['Genotype']),
        indexes=expand(["%s/{sid}_{gt}.bam.bai" % config['bcftools']['call']['idir']], zip, 
                sid = config['t']['SampleID'],
                gt = config['t']['Genotype']),
    output:
        "%s/{region,.+(:[0-9]+-[0-9]+)?}.bcf" % config['bcftools']['call']['odir']
    params:
        mpileup="--region {region}",
        call=""
    log:
        "%s/bcftools_call/{region}.log" % config['dirl']
    wrapper:
        "0.27.0/bio/bcftools/call"

