library('tidyverse')
library('cowplot')
theme_set(theme_cowplot(font_size = 7))


args <- commandArgs(trailingOnly = TRUE)
insert_dir <- args[1]
output_file <- args[2]

insert_files <- dir(insert_dir, full.names = T)

insert_df <- lapply(insert_files, function(x) {
  insert_data <- read_tsv(x, skip = 10) %>%
    dplyr::select(insert_size, count = All_Reads.fr_count) %>%
    mutate(sample_label = gsub('_insert_size_histogram.data.txt', '', str_extract(x, 'Donor_[1-4]_[[:graph:]]{0,100}')))
  total <- sum(insert_data$count)
  
  classify_fragments <- mutate(insert_data, fragment_type = case_when(
    insert_size >= 50 & insert_size <= 100 ~ 'nucleosome free',
    insert_size >= 147 & insert_size <= 2*147 ~ 'mononucleosome',
    insert_size > 2*147 ~ 'long',
    insert_size < 50 ~ 'too_short')) 
  
  fragment_summary <-  group_by(classify_fragments, sample_label, fragment_type) %>%
    summarise(fraction = sum(count)/total * 100) %>%
    filter(fragment_type %in% c('nucleosome free', 'mononucleosome'))
}) %>%
  bind_rows()

fragment_average <- group_by(insert_df, fragment_type) %>% summarise(avg = mean(fraction))

p1 <- ggplot(insert_df, aes(x = fraction, y = sample_label, fill = fragment_type)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  ylab('Sample') +
  xlab('% of fragments') +
  scale_x_continuous(n.breaks = 10) +
  ggtitle('Nucleosome Free vs. Mononucleosome Fragments') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom') +
  coord_cartesian(expand = F) +
  geom_vline(data = fragment_average, aes(xintercept = avg, color = fragment_type), lty = 'dashed')

save_plot(output_file, p1, base_height = 10, base_width = 7)
