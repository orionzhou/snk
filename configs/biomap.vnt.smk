config_default: /home/springer/zhoux379/git/snk/configs/base.yml
dir_project: /scratch.global/zhoux379/reseq/data/21_biomap
dirw: /home/springer/zhoux379/projects/reseq/data/21_biomap
dir_cache: /scratch.global/zhoux379/reseq
biomap:
  refo: /home/springer/zhoux379/data/genome/B73/download/raw.fna
  refn: /home/springer/zhoux379/data/genome/B73/21_dbs/gatk/db.fasta
  chain4: /home/springer/zhoux379/data/genome/B73/08_seq_map/mapf.chain
  ppn: 1
  runtime: 5
  mem: 20
  aruntime: 5
  amem: 10
  if01: 01.bcftools.vcf.gz
  of03: 03.liftover.vcf.gz
  of04: 04.snp.vcf.gz
  of05: 05.idl.vcf.gz
  crossmap:
    id: crsm
    ppn: 1
    runtime: 2
    mem: 30
    aruntime: 2
    amem: 10
  vcfnorm:
    id: vnorm
    ppn: 1
    runtime: 10
    mem: 30
    aruntime: 5
    amem: 10
  liftover:
    id: lif
    ppn: 1
    runtime: 5
    mem: 30
    aruntime: 5
    amem: 10
  merge:
    id: mg
    ppn: 1
    runtime: 10
    mem: 30
    aruntime: 10
    amem: 10
  filter:
    id: flt
    ppn: 1
    runtime: 15
    mem: 20
    aruntime: 5
    amem: 10
  sample:
    id: samp
    ppn: 1
    runtime: 10
    mem: 30
    aruntime: 10
    amem: 2
  phylo:
    id: phy
    q: 0
    ppn: 12
    runtime: 30
    mem: 20
    aruntime: 10
    amem: 2

