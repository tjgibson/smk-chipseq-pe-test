rule make_d4_ind:
	input:
		bam = "results/aligned_reads/filtered/{sample}.bam",
		bai = "results/aligned_reads/filtered/{sample}.bam.bai"
	output:
		"results/d4/individual/{sample}.d4"
	conda:
		"../envs/d4.yaml"
	params:
		extra=config["params"]["d4_ind"] 
	threads: 8
	shell:
		"d4tools create {params.extra} -t {threads} {input.bam} {output}"

rule merge_bam:
	input:
		get_bam_merge
	output:
		temp("results/aligned_reads/merged/{sample_group}.bam")
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
        temp("results/aligned_reads/merged/{sample}.bam.bai")
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"
        
rule make_d4_merged:
	input:
		bam = "results/aligned_reads/merged/{sample}.bam",
		bai = "results/aligned_reads/merged/{sample}.bam.bai"
	output:
		"results/d4/merged/{sample}.d4"
	conda:
		"../envs/d4.yaml"
	params:
		extra=config["params"]["d4_merged"] 
	threads: 8
	shell:
		"d4tools create {params.extra} -t {threads} {input.bam} {output}"


# rule zscore_normalize_ind_bigwigs:
# 	input:
# 		"results/bigwigs/coverage/individual/{sample}.bw"
# 	output:
# 		"results/bigwigs/zscore_normalized/individual/{sample}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/zscore_normalize_bw.R"
# 
# rule zscore_normalize_merged_bigwigs:
# 	input:
# 		"results/bigwigs/coverage/merged/{sample}.bw"
# 	output:
# 		"results/bigwigs/zscore_normalized/merged/{sample}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/zscore_normalize_bw.R"
		
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
# 		bw="results/bigwigs/zscore_normalized/individual/{sample}.bw",
# 		scaling_factors="results/scaling_factors/individual_scaling_factors.tsv"
# 	output:
# 		"results/bigwigs/spikeIn_normalized/individual/{sample}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/spikeIn_normalize_bw.R"
# 
# rule spikeIn_normalize_merged_bigwigs:
# 	input:
# 		bw="results/bigwigs/zscore_normalized/merged/{sample}.bw",
# 		scaling_factors="results/scaling_factors/merged_scaling_factors.tsv"
# 	output:
# 		"results/bigwigs/spikeIn_normalized/merged/{sample}.bw"
# 	conda:
# 		"../envs/zscore_normalize_bw.yaml"
# 	script:
# 		"../scripts/spikeIn_normalize_bw.R"