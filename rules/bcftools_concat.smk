rule bcftools_concat:
    input:
        calls=expand(["%s/{region}.bcf" % config['bcftools']['concat']['idir']], 
                region = config['bcftools']['regions'].split(" "))
    output:
        config["bcftools"]["concat"]["outfile"]
    params:
        config['bcftools']['concat']['extra']
    wrapper:
        "0.27.0/bio/bcftools/concat"
