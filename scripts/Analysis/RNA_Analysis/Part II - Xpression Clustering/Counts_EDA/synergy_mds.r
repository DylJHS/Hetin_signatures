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

# 5.  MDS
# Load Transposed transformed data 
t_scaled_log_data <- read.csv("Transposed_scaled_&_reduced_df.csv", row.names = 1)

print("t_scaled_log_data")
print(head(t_scaled_log_data[, 0:5]))

sample_dist_matrix <- dist(t_scaled_log_data)
sample_dist_df <- as.data.frame(as.matrix(sample_dist_matrix))

print("sample_dist_matrix")
print(head(sample_dist_df[, 0:5]))

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