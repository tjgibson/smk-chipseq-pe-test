# setup ------------------------------------------------------------------------
library(tidyverse)

# define input files -----------------------------------------------------------
stat_files <- snakemake@input[["mapping_stats"]]
names(stat_files) <- stat_files

# get total mapped reads for each sample  --------------------------------------
library_sizes <- stat_files |>
  map_dfr(read_tsv, .id = "source", col_names = c("chromosome", "chromosome_size", "mapped_reads", "unmapped_reads")) |>
  mutate(sample_name = gsub("_filtered.idxstats", "", basename(source))) |> 
  group_by(sample_name) |> 
  summarise(total_reads = sum(mapped_reads))

# write output table to file ---------------------------------------------------
library_sizes |> write_tsv(snakemake@output[[1]])