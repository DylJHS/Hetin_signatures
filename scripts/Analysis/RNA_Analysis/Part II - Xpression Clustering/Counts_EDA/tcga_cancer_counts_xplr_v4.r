# Load the libraries
print("Started")

library(conflicted)
library(dplyr)
library(knitr)
library(tidyverse)
library(edgeR)
library(limma)
library(data.table)
library(ggpubr)

print("Loaded packages \n")

# REUSABLE FUNCTIONS

extract_element <- function(strings, index) {
  # Split each string by "." and extract the third element
  element_list <- sapply(stringr::str_split(strings, "\\."),
                         function(x) x[index])
  return(element_list)
}

untransform <- function(x) {
  # Function to convert the log transformed counts back into original counts
  return(ceiling((2^x) - 1))
}

# Load the datasets

# RNAseq count data
original <- read.csv("tcga_gene_expected_count.csv")
test <- original[0:10000,0:15]

# Metadata
tss_meta <- read.csv("../../../data/meta/tissueSourceSite.tsv", sep = "\t")
abbrv_meta <- read.csv("../../../data/meta/bcrBatchCode.tsv", sep = "\t")
sampletype <- read.csv("../../../data/meta/sampleType.tsv",
  sep = "\t", colClasses = c("character")
)

# Gene Metadata
gene_ids <- read.delim(
  "../../../data/meta/TCGA_PanCan_TPM_Gene_Annotations.txt"
)

# SOI genes
soi <- read.csv("../../../data/TCGA_mRNA_TPM_SOI.csv")
soi_genes <-soi[,2]

print("Loaded datasets \n")


transformed_data <- test
rownames(transformed_data) <- NULL

real_count_data <- transformed_data %>%
  mutate_at(vars(-1), untransform)

# Format the data

count_data_ordered <- real_count_data[, order(colnames(real_count_data))]

# Combine the metadata
meta <- dplyr::left_join(tss_meta %>%
                           dplyr::select(c("TSS.Code", "Study.Name")) %>%
                           distinct() %>%
                           sapply(trimws) %>%
                           as.data.frame(),
                         abbrv_meta %>%
                           dplyr::select(c("Study.Abbreviation", "Study.Name")) %>% 
                           distinct()%>% 
                           sapply(trimws) %>% 
                           as.data.frame(), 
                         by = "Study.Name")

# CONVERT THE GENE IDS INTO GENE NAMES

counts_data <- dplyr::right_join(gene_ids %>% 
                    dplyr::select(c("id", "gene")) %>% 
                    sapply(trimws) %>% 
                    as.data.frame(),
                  count_data_ordered,
                  by = c("id" = "sample")) %>% 
  dplyr::select(-"id")

count_data <- counts_data %>%
  mutate(gene = trimws(gene))

print("\n count_data")
print(head(count_data[, 1:10]))


# 2.  Edit the samples (column names) into ids so as to be able to identify the participant and the sample and cancer type

# 2.1. extract the different elements from the sample names

ids <- colnames(count_data)[-1]
participants <- extract_element(ids, 3)
condition <- extract_element(ids, 4)
tissue_type <- extract_element(ids, 2)

# 2.2. Reform the column names using the extracted elements and reset the index as the Genes

column_names <- c("Genes" ,  paste(tissue_type, participants, condition, sep = "."))
colnames(count_data) <- column_names

# Convert data to a data.table for faster processing during the grouping of the duplicate genes
count_data <- data.table::setDT(count_data)

# Combine duplicate genes together using the median of the expression
grouped <- count_data[, lapply(.SD, function(x) if (length(x) > 1) ceiling(median(x, na.rm = TRUE)) else x), by = Genes, .SDcols = -"Genes"]

data <- grouped %>%
  column_to_rownames(var = "Genes")

print("\n grouped")
print(head(grouped[, 1:10]))


# Remove the unwanted sample types from the samples based on the code
codes_to_use <- c("01", "02", "03", "04", "05", "08", "09")
samples_to_use <- data %>%  dplyr::select(ends_with(codes_to_use))

# 2.3. Assign the TSS abbreviation to the column name in place of the TSS code

column_names <- colnames(samples_to_use)
prefixes <- substr(column_names,1,2)
abbrvs <- meta$Study.Abbreviation[match(prefixes, meta$TSS.Code)]
column_names_prefixes <- paste0(abbrvs, substr(column_names, 3, nchar(column_names)))
colnames(samples_to_use) <- column_names_prefixes

# 2.4. Assign the sample type abbreviation to the column name in place of the sample type code

new_column_names <- colnames(samples_to_use)
suffixes <- substr(new_column_names, nchar(new_column_names) - 1, nchar(new_column_names))
s_abbrvs <- sampletype$Short.Letter.Code[match(suffixes, sampletype$Code)]
column_names_suffxes <- paste0(substr(new_column_names, 1, nchar(new_column_names) - 2), s_abbrvs)
colnames(samples_to_use) <- column_names_suffxes

print("\n samples_to_use")
print(head(samples_to_use[, 1:10]))

# Remove 0 genes and samples
unzero_d <- samples_to_use[rowSums(samples_to_use[]) > 0, ] # samples
cols_zeros <- which(apply(unzero_d, 2, function(x) all(x == 0)))
if (length(cols_zeros) > 0) {
  reduced <- unzero_d[, -cols_zeros]
} else {
  reduced <- unzero_d
}

# convert df to matrix 
mapped_matrix <- as.matrix(reduced)

# Calculate the normalisation factor 
d0 <- edgeR::calcNormFactors(mapped_matrix, method = "TMM")

# Scale Gene counts
# Divide the samples by their corresponding normalisation factors
scaled_data <- base::sweep(samples_to_use, 2, d0, "/") # scalling the genes (2 = columns)
scaled_data <- ceiling(scaled_data) # Round up

# Filter out non SOI genes
filtered_scaled <- scaled_data[rownames(scaled_data) %in% soi_genes, ]
print("filtered data")
print(dim(filtered_scaled))

# Log Transform the counts
log_data <- apply(filtered_scaled, MARGIN = c(1,2), FUN = function(x) log(x + 1))

# 5.  MDS
# Transpose the data 
t_scaled_log_data <- t(log_data)
write.csv(t_scaled_log_data, "Transposed_scaled_&_reduced_df.csv", row.names = TRUE)

print("\n t_scaled_log_data")
print(head(t_scaled_log_data[, 1:3]))

sample_dist_matrix <- dist(t_scaled_log_data)
sample_dist_df <- as.data.frame(as.matrix(sample_dist_matrix))

print("\n sample_dist_matrix")
print(head(sample_dist_df[, 1:10]))

write.csv(sample_dist_df, "sample_dist_matrix.csv", row.names = TRUE)

# Produce the MDS object.
mds <- cmdscale(sample_dist_matrix)

# Convert the MDS object to a dataframe with additional features
mds_df <- as.data.frame(mds) %>% 
  dplyr::mutate(cancer_type = extract_element(rownames(t_scaled_log_data), 1),
         sample_type = extract_element(rownames(t_scaled_log_data), 3),
         synergy_type = paste(cancer_type, sample_type, sep = "-"))


write.csv(mds_df, "sample_mds.csv", row.names = TRUE)

# Create the plots of the MDS.
sample_mds <- ggpubr::ggscatter(mds_df, x = "V1", y = "V2",
  color = "sample_type", # Colour based on the sample type
  size = 1,
  repel = TRUE
)

cancer_mds <- ggpubr::ggscatter(mds_df, x = "V1", y = "V2", 
  color = "cancer_type", # Colour based on the cancer type
  size = 1,
  repel = TRUE) + theme(legend.position = "none")

synergy_mds <- ggpubr::ggscatter(mds_df, x = "V1", y = "V2", 
          color = "synergy_type", # Colour based on the synergy between the cancer and sample type
          size = 1,
          repel = TRUE)+
  theme(legend.position = "none")

# Save the plots.
ggsave("sample_mds_plot.png",   
       plot = sample_mds,       
       width = 8,              
       height = 6,             
       dpi = 300) 

ggsave("cancer_mds_plot.png",   
       plot = cancer_mds,       
       width = 8,              
       height = 6,             
       dpi = 300) 

ggsave("synergy_mds_plot.png",   
       plot = synergy_mds,       
       width = 8,              
       height = 6,             
       dpi = 300) 