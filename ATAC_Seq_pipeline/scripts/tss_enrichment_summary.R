library("tidyverse")
library("cowplot")
theme_set(theme_cowplot(font_size = 7))

args <- commandArgs(trailingOnly = TRUE)
enrichment_dir <- args[1]
output_file <- args[2]

enrichment_files <- dir(enrichment_dir, full.names = T, pattern = ".csv")

enrichment_df <- lapply(enrichment_files, function(x) {
  read_csv(x)
  }) %>%
  bind_rows()

p1 <- ggplot(enrichment_df, aes(x = max_tss_enrichment, y = Sample)) +
  geom_bar(stat = "identity") +
  xlab('TSS Enrichment') +
  scale_x_continuous(n.breaks = 20) +
  ggtitle('TSS Enrichment') +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(expand = F)

save_plot(output_file, p1, base_height = 10, base_width = 7)