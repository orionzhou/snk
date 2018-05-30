rule sambamba_sort:
    input:
        "%s/{pre}.bam" % config['sambamba_sort']['idir']
    output: 
        "%s/{pre}.bam" % config['sambamba_sort']['odir']
    params:
        config['sambamba']['sort']['extra']
    threads:
        config['sambamba']['threads']
    wrapper:
        "0.23.1/bio/sambamba/sort"
        
