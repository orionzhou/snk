rule star:
    input:
        fq1 = "14.trim/{sid}.fq.gz"
    output:
        "21.bam.raw/{sid}_{gt}/Aligned.out.bam"
    log:
        "logs/star/{sid}_{gt}.log"
    params:
        index = config["star"]["index"],
        outprefix = "21.bam.raw/{sid}_{gt}/",
        extra = lambda wildcards: "--outFilterType BySJout "
            "--outFilterMultimapNmax 20 "
            "--alignSJoverhangMin 8 "
            "--alignSJDBoverhangMin 1 "
            "--outFilterMismatchNmax 999 "
            "--outFilterMismatchNoverReadLmax 1.0 "
            "--alignIntronMin 20 "
            "--alignIntronMax 1000000 "
            "--alignMatesGapMax 1000000 "
            "--varVCFfile " + config['vcf'][wildcards.gt] + " "
            "--waspOutputMode SAMtag "
            #"--outSAMattributes All vA vG vW vR "
            "--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch vG vA vW "
            "--outSAMattrRGline ID:{sid} SM:{sid} "
            "--outSAMunmapped Within KeepPairs "
    threads:
        config["star"]["threads"]
    wrapper:
        "0.23.1/bio/star/align"

rule sambamba_sort:
    input:
        "21.bam.raw/{pre}/Aligned.out.bam"
    output: 
        "22.bam/{pre}.bam"
    params:
        config['sambamba_sort']['extra']
    threads:
        config['sambamba_sort']['threads']
    wrapper:
        "0.23.1/bio/sambamba/sort"

