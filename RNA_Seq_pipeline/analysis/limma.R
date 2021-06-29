library('limma')
library('edgeR')
library('tidyverse')
library('stringr')
library('cowplot')
library('viridis')
theme_set(theme_cowplot(font_size = 10))

# Read in gene counts
counts <- read_delim('~/Google_Drive/Research/analysis_data/tf_screens/rna_seq/il2ra_hits/counts/dedup_counts.txt',
                     delim = '\t', skip = 1)
colnames(counts) <- gsub('output/bam/dedup/',
                         '', gsub('.dedup.bam',
                                  '', colnames(counts)))

counts <- as.data.frame(counts)
rownames(counts) <- counts$Geneid

# Run limma for each KO vs AAVS1 controls
run_limma <- function(gene_ko) {
  # Select relevant samples and generate count matrix
  count_mat <- as.matrix(dplyr::select(counts, matches(paste0('AAVS1|', gene_ko)), -Donor_4_AAVS1_6))
  
  # Filter low count reads
  min_cpm <- ceiling(10/(min(colSums(count_mat))/1e6))
  count_mat_filtered <- count_mat[rowSums(cpm(count_mat) > min_cpm) >= 3, ]
  print(min_cpm)
  print(dim(count_mat_filtered))
  
  # Pull out just sample names and set up design matrix
  layout <- data.frame(sample = colnames(count_mat_filtered)) %>%
    mutate(ko = gsub('Donor_[0-9]_', '', sample),
           ko = gsub('AAVS1_[0-9]', 'AAVS1', ko),
           donor = str_extract(sample, 'Donor_[0-9]'))
  print(layout)
  design <- model.matrix(~ 0 + layout$ko + layout$donor)
  colnames(design) <- gsub('layout[$]ko|layout[$]donor', '', colnames(design))
  
  # Get gene names
  geneID_geneSymbol <- read_tsv('~/Google_Drive/Research/genomes/ensembl_GeneID_GeneSymbol.txt') %>%
    rename(ens_id = `Gene stable ID`, gene_name = `Gene name`)
  gene_annotation <- data.frame(ens_id = sapply(rownames(count_mat_filtered), function(y) str_split(y, '[.]')[[1]][[1]])) %>%
    left_join(., geneID_geneSymbol, by = 'ens_id')
  
  # Normalize to log2(CPM + 0.5) with voom
  d <- DGEList(count_mat_filtered, genes = gene_annotation)
  d <- calcNormFactors(d)
  v <- voom(d, design, plot = T)
  
  # Differential expression
  cmd <- paste0("contrast.matrix <- makeContrasts(", gene_ko, " - AAVS1, levels = design)")
  eval(parse(text = cmd))
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  x_results <- topTable(fit2, coef = 1, number = Inf, confint = T) %>%
    mutate(sample = gene_ko)
}

samples <- unique(str_subset(gsub('Donor_[0-9]_', '', str_subset(colnames(counts), 'Donor')), 'AAVS1', negate = TRUE))

limma_results <- lapply(samples, function(x) run_limma(x))

write_rds(limma_results, paste0('~/Google_Drive/Research/analysis_data/tf_screens/rna_seq/il2ra_hits/analysis_data/limma-', Sys.Date(), '.RDS'))

