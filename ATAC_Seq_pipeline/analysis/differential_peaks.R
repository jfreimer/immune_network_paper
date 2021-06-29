library('tidyverse')
library("vroom")
library("DESeq2")
library('edgeR')

# ------------------------------------------------------------------------------------------------
# TSS enrichment scores as surrogate for sample quality
enrichment_files <- dir("~/Google_Drive/Research/analysis_data/tf_screens/atac/il2ra_hits/tss", full.names = T, pattern = ".csv")

enrichment_df <- lapply(enrichment_files, function(x) {
  read_csv(x)
}) %>%
  bind_rows() %>%
  mutate(Sample = gsub("_[A-Z][0-9]{0,2}$", "", Sample))

# Function to run DESeq2 on each sample
run_deseq2_tss_unscaled <- function(gene_ko, count_df) {
  
  # Select relevant samples and generate count matrix
  count_mat <- as.matrix(dplyr::select(count_df, matches(paste0("AAVS1|", gene_ko))))
  
  # Filter low count reads
  min_cpm <- ceiling(10 / (min(colSums(count_mat)) / 1e6))
  count_mat_filtered <- count_mat[rowSums(cpm(count_mat) > min_cpm) >= 3, ]
  print(dim(count_mat_filtered))
  
  # Pull out just sample names and set up design matrix
  layout <- data.frame(sample = colnames(count_mat_filtered)) %>%
    mutate(
      ko = gsub("Donor_[0-9]_", "", sample),
      ko = gsub("AAVS1_[0-9]", "AAVS1", ko),
      donor = str_extract(sample, "Donor_[0-9]")
    ) %>%
    inner_join(., enrichment_df, by = c("sample" = "Sample"))
  
  print(layout)
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = count_mat_filtered,
    colData = layout,
    design = ~ max_tss_enrichment + donor + ko
  )
  
  dds <- DESeq(dds)
  res <- results(dds, c("ko", gene_ko, "AAVS1"))
  resDF <- as.data.frame(res) %>%
    mutate(
      sample = gene_ko,
      peakName = rownames(res)
    )
}


# read in counts
counts <- vroom('~/Google_Drive/Research/analysis_data/tf_screens/atac/il2ra_hits/counts/count_mat_peaks_cluster150bp_peak_size_350bp.txt')

counts <- as.data.frame(counts)

rownames(counts) <- counts$name

# remove well position
colnames(counts) <- gsub("_[A-Z][0-9]{0,2}$", "", colnames(counts))

# Extract samples from count dataframe
samples <- unique(str_subset(gsub("Donor_[0-9]_", "", str_subset(colnames(counts), "Donor")), "AAVS1", negate = TRUE))

# Run deseq2 
deseq2_results <- lapply(samples, function(x) run_deseq2_tss_unscaled(x, counts))

write_rds(deseq2_results, paste0('~/Google_Drive/Research/analysis_data/tf_screens/atac/il2ra_hits/analysis_data/differential_peaks_deseq2_', Sys.Date(), '.RDS'))
