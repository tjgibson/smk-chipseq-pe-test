# setup ------------------------------------------------------------------------
library(tidyverse)
library(ShortRead)

# define input files -----------------------------------------------------------
barcode_table_fn <- snakemake@input[["barcode_table"]]
fq_files <- snakemake@input[["fq_files"]]
sample_name <- snakemake@wildcards[["sample"]]

# read in nucleosome barcodes --------------------------------------------------
nuc_barcodes <- read_csv(barcode_table_fn)

# generate empty table of barcode counts ---------------------------------------
units <- read_tsv(unit_table)

col_names <- units$unit_name

nuc_barcodes[,col_names] <- 0

# read in raw read sequences ---------------------------------------------------
chunk_size <- 1e6

for (i in seq(col_names)) {
  fq <- input_files[i]
  
  strm <- FastqStreamer(fq, readerBlockSize = chunk_size)
  on.exit(close(stream))
  
  message(paste("starting barcode search for file:", fq))  
  
  record_count <- 0
  repeat {
    record_count <- record_count + chunk_size
    message(paste("reads processed:", record_count))
    chunk <- yield(strm)
    if (length(chunk) == 0)
      break
    ## process chunk
    for (j in seq(nrow(nuc_barcodes))) {
      chunk_count <- vcountPattern(pull(nuc_barcodes[j,"sequence"]), sread(chunk)) |> sum()
      nuc_barcodes[j,col_names[i]] <- nuc_barcodes[j,col_names[i]] + chunk_count
      
    }
  }
}

nuc_barcodes |> write_tsv(snakemake@output[[1]])

# # compute spike-in percentages and scaling factors -----------------------------
# # get total read count for each sample
# library_sizes <- tibble(file_name = input_files, total_reads = 0)
# for (i in 1:nrow(library_sizes)) {
#   file_name <- library_sizes[i, "file_name"] |> pull()
#   message(paste("computing total reads for file:", file_name))
#   library_sizes[i,"total_reads"] <-  countFastq(file_name)$records
# }
# 
# 
# # read in sample names
# units <- read_tsv("config/CR_units.tsv")
# 
# sample_names <- units |> 
#   select(unit_name, fq1, fq2) |> 
#   pivot_longer(fq1:fq2, names_to = "read_pair", values_to = "file_name") |> 
#   mutate(file_name = basename(file_name)) |> 
#   select(unit_name, file_name)
# 
# # add sample names to library sizes
# library_sizes <- library_sizes |> 
#   mutate(file_name = basename(file_name)) |> 
#   left_join(sample_names)
# 
# # coompute percent spike-in based on total counts of all barcodes for each sample
# barcodes_percentages <- nuc_barcodes |> 
#   pivot_longer(4:ncol(nuc_barcodes), names_to = "file_name", values_to = "count") |> 
#   left_join(sample_names)
# 
# scaling_factors <- barcodes_percentages |> 
#   group_by(unit_name) |> 
#   summarise(n_spike_in_reads = sum(count)) |> 
#   left_join(library_sizes) |> 
#   select(-file_name) |> 
#   distinct() |> 
#   mutate(percent_spike_in_reads = n_spike_in_reads / total_reads * 100) |> 
#   mutate(scaling_factor = 1 / percent_spike_in_reads)
# 
# 
# scaling_factors |> 
#   write_tsv("CUTandRUN/results/scaling_factors.txt")
# # compute antibody specificity from spike-in percentages -----------------------
# # compute percentage of total for each barcoded target
# AB_specificity <- barcodes_percentages |> 
#   group_by(unit_name, target) |> 
#   summarise(barcode_count = sum(count))
# 
# barcode_totals <- barcodes_percentages |> 
#   group_by(unit_name) |> 
#   summarise(total_count = sum(count))
# 
# AB_specificity <- AB_specificity |> 
#   left_join(barcode_totals) |> 
#   mutate(barcode_percent = barcode_count / total_count * 100)
# 
# # generate heatmap of specificity
# AB_specificity |> 
#   ggplot(aes(x = target, y = unit_name, fill = barcode_percent)) + geom_tile() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# # estimate CR efficiency
# AB_specificity |> 
#   select(unit_name, target, barcode_percent) |> 
#   filter(target == "H3K27me3" | target == "Unmodified") |> 
#   pivot_wider(names_from = target, values_from = barcode_percent) |> 
#   mutate(efficiency = H3K27me3 / Unmodified) |> 
#   
#   ggplot(aes(x = unit_name, y = efficiency)) + geom_bar(stat = "identity") + coord_flip()

