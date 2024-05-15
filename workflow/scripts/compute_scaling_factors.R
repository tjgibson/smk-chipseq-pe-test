# setup ------------------------------------------------------------------------
# load required packages
library(tidyverse)

# import data ------------------------------------------------------------------
# read in all relevant stat files
stat_files <- snakemake@input[["mapping_stats"]]
names(stat_files) <- stat_files

mapping_stats <- stat_files |>
  map_dfr(read_tsv, .id = "source", col_names = c("chromosome", "chromosome_size", "mapped_reads", "unmapped_reads")) |>
  mutate(sample_name = gsub("_unireads.idxstats", "", basename(source)))


# add column marking reference or spike-in chromosomes
mapping_stats <- mapping_stats |>
  mutate(reference_genome = case_when(str_detect(chromosome, "spikeIn")  ~ "spikeIn",
                                      !str_detect(chromosome, "spikeIn")  ~ "reference")) 

# compute % reads mapping to reference and spike-in genomes
mapping_stats <- mapping_stats |> 
  group_by(sample_name, reference_genome) |> 
  summarise(total_reads = sum(mapped_reads)) |>
  pivot_wider(names_from = reference_genome, values_from = total_reads)



# read in unit table
sample_table <- snakemake@config[["units"]] |> 
  read_tsv() |> 
  
  # remove redundant technical replicates from sample_table
  distinct(sample_name, .keep_all = TRUE) |> 
  select(sample_name, sample_group)


# join sample_table and mapping stats
scaling_factors <- sample_table |>
  left_join(mapping_stats, by = "sample_name")



# compute normalization factors for individual samples -------------------------
individual_scaling_factors <- scaling_factors |> 
  mutate(total_mapped_reads = reference + spikeIn) |>
  mutate(percent_reference = reference / total_mapped_reads * 100, percent_spikeIn = spikeIn / total_mapped_reads * 100) |> 
  mutate(scaling_factor = 1 / percent_spikeIn) |> 
  mutate(norm_scaling_factor = scaling_factor / max(scaling_factor))

# compute normalization factors for merged replicate samples -------------------
merged_scaling_factors <- scaling_factors |> 
  group_by(sample_group) |> 
  summarise(n_reference_reads = sum(reference), n_spikeIn_reads = sum(spikeIn)) |> 
  mutate(total_reads = n_reference_reads + n_spikeIn_reads) |>
  mutate(percent_reference = n_reference_reads / total_reads * 100, percent_spikeIn = n_spikeIn_reads / total_reads * 100) |> 
  mutate(scaling_factor = 1 / percent_spikeIn) |> 
  mutate(norm_scaling_factor = scaling_factor / max(scaling_factor)) |> 
  rename(sample_name = sample_group)

# write output files -----------------------------------------------------------
individual_scaling_factors |> 
  write_tsv(snakemake@output[[1]])

merged_scaling_factors |> 
  write_tsv(snakemake@output[[2]])
