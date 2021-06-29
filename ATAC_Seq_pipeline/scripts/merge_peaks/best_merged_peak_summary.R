library('tidyverse')
library('cowplot')
theme_set(theme_cowplot())

# -------------------------------------------------------------------------------------------------
# Read in files
args <- commandArgs(trailingOnly = TRUE)

peak_combo_counts <- str_split(args[1], ',')[[1]]

#---------------------------------------------------------------------------------

count_summary <- lapply(peak_combo_counts, function(current_file) {
  read_tsv(current_file)
}) %>%
  bind_rows()

count_summary_long <- pivot_longer(count_summary, cols = starts_with('count'))

best_combo_plot <- ggplot(count_summary_long, aes(x = factor(peak_size), y = value*100, fill = factor(aggregate_distance))) +
  geom_boxplot() +
  facet_grid(rows = vars(name), scales = 'free') +
  background_grid() +
  ylab('% of reads in peaks') +
  xlab('Peak size') +
  scale_fill_discrete(name = 'Cluster aggregation\ndistance') +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  ggtitle('Find best cluster aggregation distance and peak size')
save_plot('output/merged_peaks/plots/best_combo.pdf', best_combo_plot, base_width = 18, base_height = 12)

# Write the best 5 combos to a tsv file
best_combos <- count_summary_long %>%
  filter(name == "count_1_peak") %>%
  group_by(peak_size, aggregate_distance) %>%
  summarise(median_reads_in_peak = median(value),
            mean_reads_in_peak = mean(value)) %>%
  ungroup() %>%
  top_n(., n = 10, wt = median_reads_in_peak)

write_tsv(best_combos, path = 'output/merged_peaks/best_combo.tsv')
