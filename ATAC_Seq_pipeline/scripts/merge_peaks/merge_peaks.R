library('tidyverse')
library('cowplot')
library('GenomicRanges')
library('rtracklayer')
theme_set(theme_cowplot())

# -------------------------------------------------------------------------------------------------
# Read in files
args <- commandArgs(trailingOnly = TRUE)
summit_files <- str_split(args[1], ',')[[1]]
pileup_files <- str_split(args[2], ',')[[1]]
plot_dir <- args[3]
bed_dir <- args[4]


# summit_files <- list.files('~/Desktop/Data prototyping/peak_merging/summits', pattern = '.bed$', full.names = T)
# pileup_files <- list.files('~/Desktop/Data prototyping/peak_merging/pileup', patter = '.xls', full.names = T)
# plot_dir <- '~/Desktop/output/plots/'
# bed_dir <- '~/Desktop/output/beds/'

# ----------------------------------------------------------------------------------------------------
# Calculate reproducible summits between donors
# Count peaks as reproducible if they are less than 75 basepairs apart and have at least 10 reads covering 
# the peak summit

find_reproducible_summits <- function(sampleName, distanceCutoff = 75, minPileup = 10) {
  
  # Read in summit files for all donors for given sample
  sample_summit_files <- grep(sampleName, summit_files, value = T)
  
  individual_summits <- lapply(sample_summit_files, function(x) {
    import.bed(x) %>% keepStandardChromosomes(., pruning.mode = 'coarse')
  })
  
  # Combine summits across donors
  combined_summits <- Reduce(c, individual_summits)
  combined_summits$donor <- str_extract(string = combined_summits$name, pattern = 'Donor_[0-9]')
  
  # Read in peak summary with summit pileup information from macs
  sample_pileup_files <- grep(sampleName, pileup_files, value = T)
  
  pileup <- lapply(sample_pileup_files, function(x) {
    read_tsv(x, skip = 25) %>% filter(pileup > minPileup)
  }) %>% bind_rows()
  # Filter summits based on minimum pileup
  combined_summits <- combined_summits[combined_summits$name %in% pileup$name]
  
  
  # Cluster summits between donors to calculate reproducible summits
  # First calculate gaps between summits
  # Then take all of the gaps between gaps (i.e. clusters where peak summits are)
  sample_summit_gaps <- gaps(combined_summits)[width(gaps(combined_summits)) > distanceCutoff]
  clusters <- gaps(sample_summit_gaps)
  
  # Overlap all peak summits with clusters 
  summit_cluster_overlap <- findOverlaps(query = combined_summits, subject = clusters)
  # Identify clusters that contain a peak summit from at least 2 donors
  reproducible_cluster_ids <- tibble(cluster = subjectHits(summit_cluster_overlap),
                                     donor = combined_summits$donor[queryHits(summit_cluster_overlap)]) %>%
    group_by(cluster) %>%
    summarise(n = n_distinct(donor)) %>%
    filter(n >= 2) %>%
    select(cluster) %>%
    unlist()
  reproducible_cluster <- clusters[reproducible_cluster_ids]
  reproducible_summits <- combined_summits[queryHits(findOverlaps(query = combined_summits, subject = reproducible_cluster))]
  reproducible_summits$sample <- sampleName
  
  # Calculate distance between nearest summits and make plot
  dist_df <- tibble(distance = abs(start(combined_summits) - start(combined_summits[nearest(combined_summits)])),
                    peak1_donor = combined_summits$donor,
                    peak2_donor = combined_summits[nearest(combined_summits)]$donor) %>%
    mutate(comparison = gsub('_', '', paste(peak1_donor, 'vs', peak2_donor)))
  
  dist_plot <- ggplot(dist_df, aes(x = distance, color = comparison)) +
    geom_vline(xintercept = distanceCutoff, lty = 'dashed') +
    stat_ecdf() +
    scale_y_continuous(breaks = seq(0, 1, .1), labels = seq(0, 100, 10)) +
    scale_x_continuous(breaks = seq(0, 500, 50)) +
    coord_cartesian(xlim = c(0, 500)) +
    background_grid() +
    xlab('Distance to nearest peak summit (bp)') +
    ylab('% of peak summits in each group') +
    ggtitle(paste(sampleName, ': Distance between peak summit and its\nnearest peak summit across donors')) +
    scale_color_discrete(name = 'Peak1 Sample vs\nPeak2\n Sample:') +
    annotate("label", fill = 'white', x = 250, y = 1.05,
             label = paste('% of all summits that are reproducible between donors after clustering:',
                           round(mean(countOverlaps(combined_summits, reproducible_summits) > 0)*100, 2)))
  
  #print(dist_plot)
  save_plot(paste0(plot_dir, sampleName, '_summit_distance.pdf'), dist_plot, base_height = 10)
  return(reproducible_summits)
}

# Extract sample names from summit files
samples <- unique(gsub('_[[:alnum:]]{2,3}_summits.bed', '', gsub('Donor_[1-4]_', '', str_extract(summit_files, 'Donor_[1-4]_[[:graph:]]{0,100}'))))

# Find all reproducible summits between donors and combine
all_summits <- lapply(samples, function(x) find_reproducible_summits(x)) %>% Reduce(c, .)

#############################################################
# Output summits around IL2RA for testing
subset_range <- GRanges(seqnames = 'chr10', IRanges(5902402, 6209641))
export(subsetByOverlaps(all_summits, subset_range), paste0(bed_dir, 'IL2RA_summits.bed'))
#############################################################


# ------------------------------------------------------------------------------------------------
# Want to define peak regions for downstream analysis (counting, EDA, motif searching, etc) based on summit regions.
# First need to cluster summits together by combing all summits within some distance of each other 
# into a group. Initial clusters tend to be small regions and may have other clusters of summits nearby,
# if use all of these clusters to define peaks, will end up with a number of overlapping peaks making
# it hard to uniquely assign counts to a peak. Therefore explore aggregating clusters within different
# distances into larger combined clusters. Combine this with exploring expanding each cluster into
# different size peaks. Find the combination of cluster aggregation distance and peak size that
# maximizes the number of reads falling in unique clusters.

# -----------------------------
# Initial summit clustering

# Calculate distance between all reproducible summits from all samples merged together
all_summits_gap_distance_plot <- tibble(x = width(gaps(all_summits))) %>%
  ggplot(aes(x)) + stat_ecdf() +
  scale_y_continuous(breaks = seq(0, 1, .1), labels = seq(0, 100, 10)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(xlim = c(0, 100)) +
  background_grid() +
  xlab('Distance to nearest peak summit (bp)') +
  ylab('% of peak summits') +
  geom_vline(xintercept = 20, lty = 'dashed') +
  ggtitle('Distance beween all reproducible summits from all\nsamples merged together to define initial cluster size')
save_plot(paste0(plot_dir, 'all_summits_gap_distance.pdf'), all_summits_gap_distance_plot, base_height = 10)

# First calculate gaps between summits and only keep gaps greater than cutoff (20 bp)
# Then take all of the gaps between gaps (i.e. clusters where peak summits are)
initial_summits_gaps <- gaps(all_summits)[width(gaps(all_summits)) > 20]
initial_summit_clusters <- gaps(initial_summits_gaps)

# End up with some small clusters composed of a few outlier summits. To remove these clusters, only keep clusters 
# that contain a peak summit from at least 2 donors of the same KO sample
initial_summit_cluster_overlap <- findOverlaps(query = all_summits, subject = initial_summit_clusters)

filter_initial_cluster_ids <- tibble(cluster = subjectHits(initial_summit_cluster_overlap),
                                   donor = all_summits$donor[queryHits(initial_summit_cluster_overlap)],
                             sample = all_summits$sample[queryHits(initial_summit_cluster_overlap)]) %>%
  group_by(cluster, sample) %>%
  summarise(unique_donor_count = n_distinct(donor)) %>%
  filter(unique_donor_count >= 2) %>%
  select(cluster) %>%
  unlist() %>%
  unique()

# Initial set of clusters of peak summits
initial_summit_clusters_filtered <- initial_summit_clusters[filter_initial_cluster_ids]

#
export(initial_summit_clusters_filtered, paste0(bed_dir, 'initial_filtered_clusters.bed'))
#

# -----------------------------------------------------
# Aggregate nearby clusters
# summit_clusters_distance <- tibble(dist = mcols(distanceToNearest(resize(initial_summit_clusters_filtered, width = 1, fix = 'center')))$distance)
# ggplot(summit_clusters_distance, aes(dist)) + stat_ecdf() + coord_cartesian(xlim = c(0, 500))

# Function to aggregate clusters within cutoff_distance and then calculate average summit location
# for all summits within aggregate clusters
aggregate_clusters <- function(clusters, cutoff_distance) {
  
  # Find gaps between all clusters
  cluster_gaps <- gaps(clusters)
  
  # Only keep gaps above some cutoff size
  cluster_gaps_filtered <- cluster_gaps[width(cluster_gaps) > cutoff_distance]
  
  # Find gaps between filtered gaps (i.e. regions containing clusters)
  aggregate_clusters <- gaps(cluster_gaps_filtered)
  
  # For each aggregate cluster, take all the summits within that cluster and average their location
  summit_cluster_overlap <- findOverlaps(query = all_summits, subject = aggregate_clusters)
  
  summit_cluster_groups <- tibble(summit_chr = as.character(seqnames(all_summits[queryHits(summit_cluster_overlap)])),
                                  summit_loc = start(all_summits[queryHits(summit_cluster_overlap)]),
                                  summit_cluster = subjectHits(summit_cluster_overlap)) %>%
    group_by(., summit_chr, summit_cluster) %>%
    summarise(summits_in_cluster = n(), avg_loc = mean(summit_loc))
  
  cluster_averaged_locations <- GRanges(seqnames = summit_cluster_groups$summit_chr,
                                     IRanges(start = summit_cluster_groups$avg_loc),
                                     n_summits = summit_cluster_groups$summits_in_cluster)
  
  export(aggregate_clusters, paste0(bed_dir, 'aggregate_clusters_', cutoff_distance, 'bp.bed'))
  export(cluster_averaged_locations, paste0(bed_dir, 'aggregate_clusters_averaged_summit_', cutoff_distance, 'bp.bed'))
  
  #print(summary(width(aggregate_clusters)))
  return(cluster_averaged_locations)
}

# Aggregate each cluster over different distances to get peak centers
peak_centers <- lapply(seq(50, 200, 25), function(cluster_distance) {
  cluster_center <- aggregate_clusters(initial_summit_clusters_filtered, cutoff_distance = cluster_distance)
  cluster_center$aggregate_cluster_distance <- cluster_distance
  return(cluster_center)
})

# Resize peak centers to get list of peaks centered on averaged summit located in aggregated clusters
all_peaks <- lapply(seq(100, 500, 50), function(peak_size) {
  lapply(peak_centers, function(x) {
    peaks <- resize(x = x, width = peak_size, fix = 'center')
    peaks$peak_size = peak_size
    export(object = peaks, con = paste0(bed_dir, 'peaks_cluster_', unique(peaks$aggregate_cluster_distance), 'bp_peak_size_', peak_size, 'bp.bed'))
    return(peaks)
  })
}) %>% unlist()

write_rds(x = all_peaks, path = paste0(bed_dir, 'all_peaks.RDS'))
