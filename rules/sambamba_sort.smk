rule sambamba_sort:
    input:
        "21.bam.raw/{pre}.bam"
    output: 
        "22.bam/{pre}.bam"
    params:
        config['sambamba']['sort']['extra']
    threads:
        config['sambamba']['threads']
    wrapper:
        "0.23.1/bio/sambamba/sort"
        
