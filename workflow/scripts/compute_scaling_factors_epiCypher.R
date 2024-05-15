# setup ------------------------------------------------------------------------
library(tidyverse)

# define input files -----------------------------------------------------------
barcode_table_fn <- snakemake@input[["barcode_table"]]
unit_table_fn <- snakemake@config[["units"]]
barcode_count_fns <- snakemake@input[["barcode_counts"]]
library_sizes_fns <- snakemake@input[["library_sizes"]]

# read in input data -----------------------------------------------------------
barcode_table <- read_csv(barcode_table_fn)
unit_table <- read_tsv(unit_table_fn)
barcode_counts <- barcode_count_fns |> 
  map(read_csv, col_names = c("sample_name", "sequence", "count")) |> 
  bind_rows() |> 
  distinct()

# add barcode and sample information to count table ----------------------------
barcode_counts <- barcode_counts |> 
  left_join(barcode_table, by = "sequence") |> 
  select(target, barcode, sequence, everything()) |> 
  left_join(select(unit_table, sample_name, sample_group), by = "sample_name")

# add total library sizes to count table ---------------------------------------
library_sizes <- library_sizes_fns |> 
  read_tsv() |> 
  distinct()


# compute scaling factors ------------------------------------------------------
#individual scaling factors
scaling_factors_ind <- barcode_counts |>
  group_by(sample_name) |>
  summarise(n_spike_in_reads = sum(count)) |>
  left_join(library_sizes) |> 
  distinct() |> 
  mutate(percent_spike_in_reads = n_spike_in_reads / total_reads * 100) |>
  mutate(scaling_factor = 1 / percent_spike_in_reads) |> 
  mutate(norm_scaling_factor = scaling_factor / max(scaling_factor))

# combined replicate scaling factors
scaling_factors_comb <- barcode_counts |>
  left_join(library_sizes) |> 
  distinct() |> 
  group_by(sample_group) |>
  summarise(n_spike_in_reads = sum(count), comb_total_reads = sum(total_reads)) |>
  mutate(percent_spike_in_reads = n_spike_in_reads / comb_total_reads * 100) |>
  mutate(scaling_factor = 1 / percent_spike_in_reads) |> 
  mutate(norm_scaling_factor = scaling_factor / max(scaling_factor)) |> 
  rename(sample_name = sample_group)

# write output_data ------------------------------------------------------------
scaling_factors_ind |> write_tsv(snakemake@output[[1]])
scaling_factors_comb |> write_tsv(snakemake@output[[2]])
barcode_counts |> write_tsv(snakemake@output[[3]])
