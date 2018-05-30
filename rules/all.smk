if 'trimmomatic' in config['modules']:
    if config['paired']:
        include: "trimmomatic_pe.smk" 
    else: 
        include: "trimmomatic_se.smk"
if 'bbduk' in config['modules']:
    if not config['paired']:
        include: "bbduk_se.smk"

if 'star' in config['modules']:
    include: "star.smk"
if 'hisat2' in config['modules']:
    include: "hisat2.smk"

if 'fastqc' in config['modules']:
    include: "fastqc.smk"

if 'featurecounts' in config['modules']:
    include: "featurecounts.smk"

if 'ase' in config['modules']:
    include: "ase.smk"

if 'multiqc' in config['modules']:
    include: "multiqc.smk"
