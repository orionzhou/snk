from snakemake import shell
import vcf
input, output, params, threads, wildcards, config = snakemake.input, snakemake.output, snakemake.params, snakemake.threads, snakemake.wildcards, snakemake.config

yid, gt = wildcards.yid, wildcards.gt
adic = config['y'][yid]['ase'][gt]
inbred,sid,sid1,sid2 = adic['inbred'], adic['sid'], adic['sid1'], adic['sid2']
min_gq = params.min_gq
ref = params.ref

sids = []
if inbred:
    if gt == 'B73': sid = 'dn18a#Mo17'
    sids.append(sid)
else:
    if sid1.endswith('B73'):
        sids.append(sid2)
    elif sid2.endswith('B73'):
        sids.append(sid1)
    else:
        sids += [sid1, sid2]

sid_str = ",".join(sids)
if len(sids) == 1:
    shell("""
        bcftools view -a -s {sid_str} {input.vcf} -Ou | \
            bcftools view -i 'ALT!="*" && GQ>={min_gq} && GT="AA"' -Ou | \
            bcftools norm -f {ref} -Oz -o {output.vcf}
    """)
else:
    shell("""
        bcftools view -a -s {sid_str} {input.vcf} -Ou | \
            bcftools view -i 'ALT!="*" && N_PASS(GQ>={min_gq}) = 2 && N_PASS(GT="AA") = 1 & N_PASS(GT="RR") = 1' -Ou | \
            bcftools norm -f {ref} -Oz -o {output.vcf}
    """)

eff = output.stat.replace('.txt', '.tsv')
shell("bcftools index -t {output.vcf}")
shell("bcftools stats {output.vcf} > {output.stat}")
shell("bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' {output.vcf} | gzip -c > {eff}.gz")

vcfi = vcf.Reader(filename = output.vcf)

fho = open(params.vcf2, "w")
fho.write("##fileformat=VCFv4.2\n")
fho.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
for chrom, contig in vcfi.contigs.items():
    length = contig.length
    fho.write("##contig=<ID=%s,length=%d>\n" % (chrom, length))
fho.write("\t".join("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample".split())+"\n")

for rcd in vcfi:
    ary = [rcd.CHROM, rcd.POS, rcd.ID, rcd.REF, rcd.ALT[0], 999, '.', '.', 'GT']
    ary = [str(x) for x in ary]
    if inbred:
        call = rcd.genotype(sid)
        gty, gq = call.gt_type, call['GQ']
        assert gty == 2 and gq is not None and gq >= min_gq, 'error gt'
        fho.write("\t".join(ary + ["0|1"])+"\n")
    else:
        if sid1.endswith('B73'):
            call2 = rcd.genotype(sid2)
            gty2, gq2 = call2.gt_type, call2['GQ']
            assert gty2 == 2 and gq2 is not None and gq2 >= min_gq, 'error gt'
            fho.write("\t".join(ary + ["0|1"])+"\n")
        elif sid2.endswith('B73'):
            call1 = rcd.genotype(sid1)
            gty1, gq1 = call1.gt_type, call1['GQ']
            assert gty1 == 2 and gq1 is not None and gq1 >= min_gq, 'error gt'
            fho.write("\t".join(ary + ["1|0"])+"\n")
        else:
            call1, call2 = rcd.genotype(sid1), rcd.genotype(sid2)
            gty1, gq1 = call1.gt_type, call1['GQ']
            gty2, gq2 = call2.gt_type, call2['GQ']
            assert gq1 is not None and gq1 >= min_gq, 'error gt'
            assert gq2 is not None and gq2 >= min_gq, 'error gt'
            if gty1 == 0 and gty2 == 2:
                fho.write("\t".join(ary + ["0|1"])+"\n")
            elif gty1 == 2 and gty2 == 0:
                fho.write("\t".join(ary + ["1|0"])+"\n")

fho.close()

shell("""
    bgzip {params.vcf2}
    bcftools view -v snps -Ob {output.vcf2} -o {output.bcf}
    bcftools index -c {output.bcf}
    """)


