rule sambamba_merge:
    input:
        "21.bam.raw/{pre}.PE.bam",
        "21.bam.raw/{pre}.SE.bam",
    output:
        "22.bam/{pre}.merged.bam"
    params:
        config["sambamba"]["merge"]["extra"]
    threads:
        config["sambamba"]["threads"]
    run:
        "0.23.1/bio/sambamba/merge"
        
