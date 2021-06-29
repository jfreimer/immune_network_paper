library('GenomicRanges')
library('tidyverse')
library('cowplot')
library('zoo')
theme_set(theme_cowplot(font_size = 10))

args <- commandArgs(trailingOnly = TRUE)

tss_file <- args[1]
atac_file <- args[2]
window_size <- as.numeric(args[3])
sample <- args[4]
plot_file <- args[5]
enrichment_file <- args[6]

# tss_file <- '~/Desktop/tss_comparison/gencode_tss_hg38.bed'
# atac_file <- '~/Desktop/tss_comparison/Donor_1_AAVS1_1_C3.insertions.bed'
# window_size <- 2000
# sample <- 'test'
# plot_file <- 'test.pdf'
# enrichment_file <- 'test_enrichment.csv'

# Read in TSS bed file and make TSS GRanges
tss_bed <- read_tsv(tss_file, col_names = F)
tss <- GRanges(seqnames = tss_bed$X1, IRanges(start = tss_bed$X2, end = tss_bed$X3), strand = tss_bed$X6)
tss <- keepStandardChromosomes(tss, pruning.mode = 'coarse')

# Read in ATAC-Seq read bed file and make ATAC-Seq reads GRanges
# Filter atac_file to only use chromosomes present in bed file
atac_bed <- read_tsv(atac_file, col_names = F)
atac <- GRanges(seqnames = atac_bed$X1, IRanges(start = atac_bed$X2, end = atac_bed$X3), strand = '*')
atac <- keepStandardChromosomes(atac, pruning.mode = 'coarse')

# Find overlaps between atac reads and tss within window size
overlap <- findOverlaps(query = tss, subject = atac, maxgap = window_size, ignore.strand = TRUE)

# Get + strand and - strand TSS
tss_minus <- which(strand(tss[queryHits(overlap)]) == '-')
tss_plus <- which(strand(tss[queryHits(overlap)]) == '+')

# Calculate distance between atac read and TSS taking into account strand
distance_minus <- start(tss[queryHits(overlap[tss_minus])]) - start(atac[subjectHits(overlap[tss_minus])])
distance_plus <- start(atac[subjectHits(overlap[tss_plus])]) - start(tss[queryHits(overlap[tss_plus])])

# Calculate number of reads at each distance to the TSS within the set window
read_coverage <- tibble(dist_to_tss = c(distance_minus, distance_plus)) %>%
  group_by(dist_to_tss) %>%
  summarise(read_count = n()) %>%
  filter(dist_to_tss >= -window_size, dist_to_tss <= window_size)

# Normalize number of reads at each base by average background in first and last 100 bases of window
norm_factor <- dplyr::filter(read_coverage, dist_to_tss <= -window_size+100 | dist_to_tss >= window_size-100) %>%
  summarise(avg = mean(read_count)) %>%
  unlist()


read_coverage <- mutate(read_coverage,
                        tss_enrichment = read_count/norm_factor,
                        tss_enrichment_smooth = rollmean(tss_enrichment, 51, 1))

tss_plot <- ggplot(read_coverage) +
  geom_point(aes(x = dist_to_tss, y = tss_enrichment), color = 'lightgrey') +
  geom_line(aes(x = dist_to_tss, y = tss_enrichment_smooth), size = 1.5) +
  scale_y_continuous(n.breaks = 10) +
  ylab('TSS Enrichment') +
  xlab('Distance to TSS') +
  ggtitle(sample) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(expand = F)

save_plot(plot_file, tss_plot)

write_csv(data.frame(Sample = sample, max_tss_enrichment = max(read_coverage$tss_enrichment_smooth)), enrichment_file)
