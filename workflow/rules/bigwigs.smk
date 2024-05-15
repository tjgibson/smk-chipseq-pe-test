# rule make_bigwigs_ind:
# 	input:
#		bam = "results/aligned_reads/filtered/{sample}.bam",
#		bai = "results/aligned_reads/filtered/{sample}.bam.bai"
#	output:
#		"results/bigwigs/coverage/individual/{sample}.bw"
#	conda:
#		"../envs/deeptools.yaml"
#	params:
#		extra=config["params"]["bigwigs_ind"] 
#	threads: 8
#	shell:
#		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"

rule merge_bam:
	input:
		get_bam_merge
	output:
		"results/aligned_reads/merged/{sample_group}.bam"
	params:
		"" # optional additional parameters as string
	threads:  # Samtools takes additional threads through its option -@
		8     # This value - 1 will be sent to -@
	wrapper:
		"v1.1.0/bio/samtools/merge"		

rule samtools_index_merged:
    input:
        "results/aligned_reads/merged/{sample}.bam"
    output:
        "results/aligned_reads/merged/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"
        
#rule make_bigwigs_merged:
#	input:
#		bam = "results/aligned_reads/merged/{sample}.bam",
#		bai = "results/aligned_reads/merged/{sample}.bam.bai"
#	output:
#		"results/bigwigs/coverage/merged/{sample}.bw"
#	conda:
#		"../envs/deeptools.yaml"
#	params:
#		extra=config["params"]["bigwigs_merged"] 
#	threads: 8
#	shell:
#		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"

rule input_normalize_ind_bigwigs:
	input:
		unpack(get_bamCompare_input_ind)
	output:
		"results/bigwigs/input_normalized/individual/{sample}.bw"
	conda:
		"../envs/deeptools.yaml"
	params:
		extra=config["params"]["input_normalize_bigwigs_ind"] 
	threads: 8
	shell:
		"bamCompare -b1 {input.treatment} -b2 {input.control} -o {output} -p {threads} {params.extra}"

rule input_normalize_merged_bigwigs:
	input:
		unpack(get_bamCompare_input_merged)
	output:
		"results/bigwigs/input_normalized/merged/{sample_group}.bw"
	conda:
		"../envs/deeptools.yaml"
	params:
		extra=config["params"]["input_normalize_bigwigs_merged"] 
	threads: 8
	shell:
		"bamCompare -b1 {input.treatment} -b2 {input.control} -o {output} -p {threads} {params.extra}"

# rule zscore_normalize_ind_bigwigs:
# 	input:
# 		"results/bigwigs/coverage/individual/{sample}_{frag_size}.bw"
# 	output:
# 		"results/bigwigs/zscore_normalized/individual/{sample}_{frag_size}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/zscore_normalize_bw.R"
# 
# rule zscore_normalize_merged_bigwigs:
# 	input:
# 		"results/bigwigs/coverage/merged/{sample}_{frag_size}.bw"
# 	output:
# 		"results/bigwigs/zscore_normalized/merged/{sample}_{frag_size}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/zscore_normalize_bw.R"
# 		
# if config["use_spikeIn"] and not config["epiCypher_spikeIn"]:
# 	rule compute_scaling_factors_genome:
# 		input:
# 			mapping_stats=get_scaling_input
# 		output:
# 			"results/scaling_factors/individual_scaling_factors.tsv",
# 			"results/scaling_factors/merged_scaling_factors.tsv"
# 		conda:
# 			"../envs/zscore_normalize_bw.yaml"
# 		script:
# 			"../scripts/compute_scaling_factors.R"
# 
# 
# if config["use_spikeIn"] and config["epiCypher_spikeIn"]:
# 	rule count_epiCypher_barcodes:
# 		input:
# 			get_NGmerge_input,
# 		params:
# 			barcode="{barcode_sequence}"
# 		conda:
# 			"../envs/fqgrep.yaml"
# 		output:
# 			temp("results/scaling_factors/{sample}_{barcode_sequence}_count.csv"),
# 		shell:
# 			"""
# 			count=$(fqgrep -c {params.barcode} {input})
# 			echo {wildcards.sample},{params.barcode},$count > {output}
# 			"""
# 
# 	rule compute_library_sizes:
# 		input:
# 			mapping_stats=get_mapping_stat_fns,
# 		output:
# 			temp("results/scaling_factors/total_mapped_reads.tsv"),
# 		conda:
# 			"../envs/zscore_normalize_bw.yaml"
# 		script:
# 			"../scripts/compute_library_sizes.R"
# 
# 	rule compute_scaling_factors_epiCypher:
# 		input:
# 			barcode_table=config["epiCypher_barcodes"],
# 			units=config["units"],
# 			barcode_counts=get_scaling_input_epiCypher,
# 			library_sizes = "results/scaling_factors/total_mapped_reads.tsv",
# 		output:
# 			"results/scaling_factors/individual_scaling_factors.tsv",
# 			"results/scaling_factors/merged_scaling_factors.tsv",
# 			"results/scaling_factors/epiCypher_barcode_counts.tsv",
# 		conda:
# 			"../envs/zscore_normalize_bw.yaml"
# 		script:
# 			"../scripts/compute_scaling_factors_epiCypher.R"
# 		
# rule spikeIn_normalize_ind_bigwigs:
# 	input:
# 		bw="results/bigwigs/zscore_normalized/individual/{sample}_{frag_size}.bw",
# 		scaling_factors="results/scaling_factors/individual_scaling_factors.tsv"
# 	output:
# 		"results/bigwigs/spikeIn_normalized/individual/{sample}_{frag_size}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/spikeIn_normalize_bw.R"
# 
# rule spikeIn_normalize_merged_bigwigs:
# 	input:
# 		bw="results/bigwigs/zscore_normalized/merged/{sample}_{frag_size}.bw",
# 		scaling_factors="results/scaling_factors/merged_scaling_factors.tsv"
# 	output:
# 		"results/bigwigs/spikeIn_normalized/merged/{sample}_{frag_size}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/spikeIn_normalize_bw.R"
