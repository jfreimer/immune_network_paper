library('tidyverse')
library('GenomicRanges')

# -------------------------------------------------------------------------------------------------
# Read in files
args <- commandArgs(trailingOnly = TRUE)

all_peaks <- read_rds(args[1])
insertion_file <- args[2]

 
# all_peaks <- read_rds()/oak/stanford/groups/pritch/users/jake/tfscreens/atac/il2ra_hits/output/merged_peaks/beds
# insertion_files <-

# -----------------------------------------------------
# Count how many reads (ATAC-Seq insertion sites) overlap a single peak or multiple peaks
# for different peak sizes
read_overlap_count <- function(reads, peaks) {
  read_count <- countOverlaps(reads, peaks)
  read_count_summary <- tibble(peak_size = unique(peaks$peak_size),
                               aggregate_distance = unique(peaks$aggregate_cluster_distance),
                               count_any_peak = mean(read_count > 0),
                               count_1_peak = mean(read_count == 1),
                               count_multi_peaks = mean(read_count > 1))
}


# -----------------------------------------------------
# Summarize best aggregation and peak size combination for all samples
current_sample <- gsub('.insertions.bed', '', str_extract(insertion_file, 'Donor_[1-4]_[[:graph:]]{0,100}'))

insertions_df <- vroom::vroom(insertion_file, delim = '\t',
                              col_names = c('chr', 'start', 'stop', 'readid', 'score', 'strand'),
                              col_select = c(chr, start))
insertions <- GRanges(seqnames = insertions_df$chr, IRanges(start = insertions_df$start+1, width = 1))
rm(insertions_df)

sample_count_summary <- lapply(all_peaks, function(x) read_overlap_count(reads = insertions, peaks = x)) %>%
  bind_rows() %>%
  mutate(sample = current_sample)

write_tsv(x = sample_count_summary, path = paste0('output/merged_peaks/best_combo_count/', current_sample, '_peaks_count.tsv'))
