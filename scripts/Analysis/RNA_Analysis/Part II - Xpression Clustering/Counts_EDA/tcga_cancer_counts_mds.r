# This script is part of the counts exploraiton script and is intended to run the MDS for the data.

# Load the packages
library(plyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(conflicted)

extract_element <- function(strings, index) {
  # Split each string by "." and extract the third element
  element_list <- sapply(strsplit(strings, "\\."), function(x) x[index])
  return(element_list)
}

print("packages loaded")

# Load the transposed data.
t_data <- read.csv("Transposed_scaled_&_reduced_df.csv")
# t_data <- t_data[0:5000,0:5000]
print("t_data")
print(head(t_data[, 1:10]))
print(dim(t_data))

soi <- read.csv("../../../data/TCGA_mRNA_TPM_SOI.csv")
soi_genes <-soi[,2]

# Filter the data down to the Set of interest.
intermediate_data <- t_data %>% column_to_rownames(var = "X")
filtered_transposed <- intermediate_data[, colnames(intermediate_data) %in% soi_genes]
print("filtered data")
print(dim(filtered_transposed))

# Create the distance matrix on the filtered data.
dist_matrix <- filtered_transposed %>%
  dist() 

write.csv(dist_matrix, "sample_dist_matrix.csv", row.names = TRUE)

# Produce the MDS object.
mds <- cmdscale(dist_matrix)

# Convert the MDS object to a dataframe with additional features
mds_df <- as.data.frame(mds) %>% 
  dplyr::mutate(cancer_type = extract_element(rownames(filtered_transposed), 1),
         sample_type = extract_element(rownames(filtered_transposed), 3))


write.csv(mds_df, "sample_mds.csv", row.names = TRUE)

# Create the plots of the MDS.
sample_mds <- ggscatter(mds_df, x = "V1", y = "V2", 
          color = "sample_type", # Colour based on the sample type
          size = 1,
          repel = TRUE)

cancer_mds <- ggscatter(mds_df, x = "V1", y = "V2", 
          color = "cancer_type", # Colour based on the cancer type
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