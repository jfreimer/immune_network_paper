library('tidyverse')
library('rtracklayer')
library('GenomicRanges')
library('GenomicAlignments')

# -------------------------------------------------------------------------------------------------
# Read in files
args <- commandArgs(trailingOnly = TRUE)

final_peaks_file <- args[1]
insertion_files <- str_split(args[2], ',')[[1]]
output_file <- args[3]

# start <- Sys.time()
# final_peaks_file <- "~/Desktop/best_peaks/merged_peaks/beds/peaks_cluster_125bp_peak_size_300bp.bed"
# insertion_files <- list.files('~/Desktop/Data prototyping/insertions/', pattern = '.bed$', full.names = T)
# output_file <- "~/Desktop/count_mat_peaks_cluster125bp_peak_size_300bp.txt"
#-----------------------------------
# Import final peaks and give each peak a name
final_peaks <- import.bed(final_peaks_file)
names(final_peaks) <- paste0('peak', seq(1:length(final_peaks)))
final_peaks$name <- names(final_peaks)

# -----------------------------------------------------
# Function to count atac reads over each peak
# Won't assign atac reads overlapping more than 1 peak
get_count_matrix <- function(insertion_file) {
  insertions_df <- vroom::vroom(insertion_file, delim = '\t',
                                col_names = c('chr', 'start', 'stop', 'readid', 'score', 'strand'),
                                col_select = c(chr, start))
  insertions <- GRanges(seqnames = insertions_df$chr, IRanges(start = insertions_df$start+1, width = 1))
  rm(insertions_df)
  
  current_sample <- gsub('.insertions.bed', '', str_extract(insertion_file, 'Donor_[1-4]_[[:graph:]]{0,100}'))
  
  so <- summarizeOverlaps(final_peaks, insertions)
  count_mat <- assays(so)$counts
  colnames(count_mat) <- current_sample
  return(count_mat)
}

peak_counts <- lapply(insertion_files, get_count_matrix)

# Combine each count into 1 matrix and annotate peak information
combined_count_mat <- Reduce(cbind, peak_counts) %>% as_tibble(rownames = 'peakName')
annotated_count_mat <- final_peaks %>% as_tibble() %>% select(-score) %>% inner_join(., combined_count_mat, by = c('name' = 'peakName'))

write_tsv(annotated_count_mat, output_file)
