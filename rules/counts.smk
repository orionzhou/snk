from os import path

import numpy as np

#############
# FUNCTIONS #
#############

def normalize_counts(counts):
	"""Normalizes expression counts using DESeq's median-of-ratios approach."""

	with np.errstate(divide="ignore"):
		size_factors = estimate_size_factors(counts)
		return counts / size_factors

def estimate_size_factors(counts):
	"""Calculate size factors for DESeq's median-of-ratios normalization."""

	def _estimate_size_factors_col(counts, log_geo_means):
		log_counts = np.log(counts)
		mask = np.isfinite(log_geo_means) & (counts > 0)
		return np.exp(np.median((log_counts - log_geo_means)[mask]))

	log_geo_means = np.mean(np.log(counts), axis=1)
	size_factors = np.apply_along_axis(
		_estimate_size_factors_col, axis=0,
		arr=counts, log_geo_means=log_geo_means)

	return size_factors

#########
# RULES #
#########
def feature_counts_extra(wildcards):
	extra = config["feature_counts"]["extra"]
	if is_paired:
		extra += " -p"
	return extra

rule feature_counts:
    input:
        bam="{reference}/bam/final/{sample}.bam"
    output:
        counts="{reference}/counts/per_sample/{sample}_gene.txt",
        summary="{reference}/qc/feature_counts/{sample}_gene.txt"
    params:
        annotation=lambda wildcards: config["feature_counts"]["annotation"][wildcards.reference],
        extra=feature_counts_extra
    threads:
        config["feature_counts"]["threads"]
    log:
        "{reference}/logs/feature_counts/{sample}.txt"
    run:
        shell("featureCounts {params.extra} -a {params.annotation} -o {output.counts} -T {threads} {input.bam} > {log}")
        # Move summary to expected location
        summary_path = output.counts + ".summary"
        if summary_path != output.summary:
            shell("mv {summary_path} {output.summary}")


#def merge_inputs(wildcards):
#    all=get_samples_with_meta(wildcards.sample)
#    file_paths = ["{}/bam/sorted/{}.bam".format(wildcards.reference, samp) for samp in all]
#    print(file_paths)
#    return file_paths

def get_merge_inputs(wildcards):
    samples=get_samples()
    test_paths = ["{}/counts/per_sample/{}_gene.txt".format(wildcards.reference, samp) for samp in samples]
    return test_paths

rule merge_counts:
    input:
        get_merge_inputs
        # To merge across both references:
        #expand("{reference}/counts/per_sample/{sample}.txt", sample=get_samples(), reference=references)
    output:
        "{reference}/counts/merged.txt"
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1,index_col=list(range(6))) for fp in input)
        merged = pd.concat(frames, axis=1)

        # Extract sample names.
        merged = merged.rename(
            columns=lambda c: path.splitext(path.basename(c))[0])

        merged.to_csv(output[0], sep="\t", index=True)


rule normalize_counts:
    input:
        "{reference}/counts/merged.txt"
    output:
        "{reference}/counts/merged.log2.txt"
    run:
        counts = pd.read_csv(input[0], sep="\t", index_col=list(range(6)))
        norm_counts = np.log2(normalize_counts(counts) + 1)
        norm_counts.to_csv(output[0], sep="\t", index=True)
