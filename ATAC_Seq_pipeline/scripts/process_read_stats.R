library('tidyverse')
library('cowplot')
theme_set(theme_cowplot(font_size = 7))

# Get filenames and path from snakemake pipeline
args <- commandArgs(trailingOnly = TRUE)
read_stats_dir <- args[1]
final_reads_plot_out <- args[2]
mito_percent_plot_out <- args[3]
count_summary_plot_out <- args[4]

# Mitochondrial counts
mito_files <- dir(read_stats_dir, full.names = T, pattern = 'mito')
mito_counts <- lapply(mito_files, function(x) {
  read_csv(x)
}) %>% bind_rows()

# Mapped, unfiltered counts
mapped_files <- dir(read_stats_dir, full.names = T, pattern = 'mapped')
mapped_counts <- lapply(mapped_files, function(x) {
  read_csv(x)
}) %>% bind_rows()

# Final counts after all filtering
final_files <- dir(read_stats_dir, full.names = T, pattern = 'final')
final_counts <- lapply(final_files, function(x) {
  read_csv(x)
}) %>% bind_rows()

# Reads filtered out due to qc, 
filtered_files <- dir(read_stats_dir, full.names = T, pattern = 'filtered')
filtered_counts <- lapply(filtered_files, function(x) {
  read_csv(x)
}) %>% bind_rows()


# Final counts after all filtering
duplicated_files <- dir(read_stats_dir, full.names = T, pattern = 'duplicated')
duplicated_counts <- lapply(duplicated_files, function(x) {
  read_csv(x)
}) %>% bind_rows()

# Join count summaries into 1 dataframe
count_table <- Reduce(function(...) inner_join(..., by = 'Sample'),
                      list(mito_counts, mapped_counts, final_counts, filtered_counts, duplicated_counts))

# Plot final read count after all filtering
final_reads_plot <- ggplot(count_table, aes(x = final_reads, y = Sample)) +
  geom_bar(stat = 'identity') +
  xlab('Final Read Count') +
  scale_x_continuous(n.breaks = 10) +
  ggtitle('Final Read Count') +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(expand = F)

save_plot(final_reads_plot_out, final_reads_plot, base_height = 10, base_width = 7)

# Plot percent of mitochondrial reads in mapped libraries
mito_percent_plot <- count_table %>% mutate(mito_percent = mito_reads/mapped_reads*100) %>%
  ggplot(aes(x = mito_percent, y = Sample)) +
  geom_bar(stat = 'identity') +
  xlab('% Mitochonrial Reads') +
  scale_x_continuous(n.breaks = 20) +
  ggtitle('Mitochondrial Reads') +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(expand = F)

save_plot(mito_percent_plot_out, mito_percent_plot, base_height = 10, base_width = 7)

# Plot summary breakdown of different read % in library
count_summary_plot <- count_table %>% dplyr::select(-mapped_reads) %>%
  pivot_longer(matches('reads')) %>%
  group_by(Sample) %>%
  mutate(norm_value = value/sum(value)*100) %>%
  mutate(name_ordered = factor(name, levels = c('mito_reads', 'filtered_reads', 'duplicated_reads', 'final_reads')))  %>%
  ggplot(aes(x = norm_value, y = Sample, fill = name_ordered)) +
  geom_bar(stat = 'identity') +
  xlab('% of mapped reads') +
  ggtitle('% Reads in library') +
  scale_x_continuous(n.breaks = 20) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom') +
  scale_fill_discrete(name = 'Read Type:', labels = c("Mitochondrial", "Filtered", "Duplicated", "Final")) +
  coord_cartesian(expand = F)

save_plot(count_summary_plot_out, count_summary_plot, base_height = 10, base_width = 7)
