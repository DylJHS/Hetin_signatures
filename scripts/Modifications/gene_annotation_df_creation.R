#This script imports a GTF file containing human gene annotations, 
#extracts unique gene IDs, gene names, and chromosome names, 
# removes entries with missing values, and saves the cleaned data frame as a CSV file.

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# Load the gene annotations 
gtf_file <- "HetIN_signatures/data/secondary_data/Homo_sapiens_GRCh38_114_chr.gtf"
gtf <- import(gtf_file)

gene_df <- data.frame(
  gene_id = mcols(gtf)$gene_id,
  gene_name = mcols(gtf)$gene_name,
  seqnames = as.character(seqnames(gtf))
) %>% unique()

# Drop the NA
gene_df <- gene_df %>% drop_na(gene_id, gene_name)


# Save the df 
# write.csv(gene_df, "data/secondary_data/gene_annotations.csv", row.names = FALSE)
