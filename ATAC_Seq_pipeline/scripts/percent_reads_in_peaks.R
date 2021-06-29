library('tidyverse')
library('cowplot')
library('ggrepel')
theme_set(theme_cowplot(font_size = 7))

# Get filenames and path from snakemake pipeline
args <- commandArgs(trailingOnly = TRUE)
stats_dir <- args[1]
plot_out <- args[2]

# Percent of reads in peak
peak_combo_counts <- list.files(stats_dir, full.names = T)
count_summary <- lapply(peak_combo_counts, function(current_file) {
  read_tsv(current_file)
}) %>%
  bind_rows()

count_summary_long <- pivot_longer(count_summary, cols = starts_with('count')) %>%
  filter(peak_size == 350, aggregate_distance == 150, name == 'count_any_peak')

p <- ggplot(count_summary_long, aes(x = value*100, y = sample)) + geom_bar(stat = 'identity') +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  xlab('% of reads in peaks') +
  ggtitle('Percent of reads in peaks')

save_plot(plot_out, p, base_height = 10, base_width = 7)
