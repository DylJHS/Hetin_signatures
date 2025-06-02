library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(paletteer)


# Set the working directory. 

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship:Research Proj/Epigenetics/Janssen Group - UMC Utrecht/Main Project/Scripts")


# Load the RNA metrics dataset and remove the useless columns. 

# The data shows the samples as rows and descriptive characteristics 
# including RNA Mean, RNA Median, the Tissue Source Site (TSS) and the 
# Histology grouping as columns.

data <- fread("../Data/RNA/Full_RNA_Sample_metrics.csv")



# Create a column for the difference between the means and same for the medians . 

data[, ':=' ('Mean_diff' = Mean_set - Mean_cntrl, 'Median_diff' = Median_set - Median_cntrl)]

order_data <- setDF(data) %>% 
  arrange( Mean_diff) %>% 
  mutate(Sample = factor(Sample, levels = unique(Sample)))




# Part I. Inter-genomic expression analysis 


# Q1. Create a histogram of the differences of the average log2 expression as 
# evidence of normal distribution necessary for paired t-test.

Histo_RNA_DifMean <- ggplot(data, aes(x = Mean_diff)) + 
                          geom_histogram(colour = 1, 
                                         fill = "#49953D", 
                                         bins = 28) +
  labs(y = "Counts", x = "Differential Log2 Expression",
       title = "Counts of Differential Log2 Expression \n Across Samples") +
  theme(axis.title = element_text(size = 13,
                                  color = "black",
                                  face = "bold"),
        axis.text = element_text(color = "black",
                                 size = 10),
        panel.background = element_rect(fill = "#FFFFF0"),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        panel.border = element_rect(fill = "transparent",
                                    color = 1,
                                    linewidth = 1),
        panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = 2),
        panel.grid.major.x = element_line(color = "black", size = 0.5, linetype = 2),
        plot.title.position = "plot")
  


# Q2. Create a barplot of the samples vs their differences in the mean log2 
# RNA values coloured by TSS.

Sample_RNA_Dif_TSS<- ggplot(order_data, aes(x = Sample, y = Mean_diff, fill = TSS)) + 
  geom_col() +
  theme_minimal() +
  labs(y = "Δ Mean Log2 Expression",
       x = "Samples",
      title = "Expression Differential by Sample ")+
  scale_y_continuous(breaks = c(0, 1, 2, 3))+
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14), 
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 14), 
    legend.position = "right",
    legend.text = element_text(size= 11),
    legend.title =element_text(hjust = 0.5,
                               size = 13),
    plot.title = element_text(hjust = 0.4,
                              size = 18),
    plot.title.position = "plot",
    panel.background = element_rect(fill = "#FFFFF0"),
    panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = 2),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(size = 0.35, linetype = "solid",
                             colour = "black"),
    panel.border = element_rect(fill = "transparent",
                                color = 1,
                                linewidth = 1)) +
  guides(fill = guide_legend(ncol = 1, 
                             title = "Source Tissue Type",
                             title.position = "top")) 



# Q2. Create a barplot of the samples vs their differences in the mean 
# log2 RNA values coloured by Hist

Sample_RNA_Dif_Hist <- ggplot(order_data, aes(x = Sample, y = Mean_diff, fill = Hist_group)) + 
  geom_col() +
  theme_minimal() +
  labs(y = "Δ Mean Log2 Expression",
       x = "Samples",
       title = "Expression Differential by Sample ")+
  scale_y_continuous(breaks = c(0, 1, 2, 3))+
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14), 
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 14), 
    legend.position = "right",
    legend.text = element_text(size= 11),
    legend.title =element_text(hjust = 0.5,
                               size = 13),
    plot.title = element_text(hjust = 0.4,
                              size = 18),
    plot.title.position = "plot",
    panel.background = element_rect(fill = "#FFFFF0"),
    panel.grid.major.y = element_line(color = "black", size = 0.5, linetype = 2),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(size = 0.35, linetype = "solid",
                             colour = "black"),
    panel.border = element_rect(fill = "transparent",
                                color = 1,
                                linewidth = 1)) +
  guides(fill = guide_legend(ncol = 1, 
                             title = "Histology Grouping",
                             title.position = "top"))

# Q1 & Q2. Create the functions that are going to perform the paired 
# t-test for each cancer type on the mean and the median stats.

# Function create a df that calculates the mean difference over the provided set
# and applies a paired t-test returning mean difference as well as the p-value.

perform_mean_test <- function(data) {
  test_result <- t.test(data$Mean_set, data$Mean_cntrl, paired = TRUE)
  return(data.frame(mean_diff = mean(data$Mean_set - data$Mean_cntrl),
                    p_value = test_result$p.value))
}



# Applies the "perform_mean_test" function seperately on each TSS group.

mean_results <- data %>%
  group_by(TSS, Hist_group) %>%
  filter(n() > 1 & !is.na(Mean_set) & !is.na(Mean_cntrl)) %>%
  do(perform_mean_test(.))



# Correct the p-values due to mutliple testing using the 
# Benjamini-Hochberg method for FDR.

results_mean <- mean_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) 



# Replace p-values with a value of 0 with a very small number to avoid
# downstream errors that would occur when applying log transformations.

min_nonzero_p <- min(results_mean$p_adjusted[results_mean$p_adjusted > 0])
results_mean$p_value_adjusted <- ifelse(results_mean$p_adjusted == 0, min_nonzero_p, results_mean$p_adjusted)



# Apply a -log10 transformation to the p-values for readability on graphs.

results_mean$neg_log_pvalue <- -log10(results_mean$p_value_adjusted)



# Sort by Histology group for readability/interpretation of graph and legend.

results_mean_complete <- results_mean %>% arrange(Hist_group) 
unique_tss <- unique(results_mean_complete$TSS)
results_mean_complete$TSS <- factor(results_mean_complete$TSS, levels = unique_tss)



# Add significance threshold. 

threshold <- -log10(0.01)



# Q1 & Q2. Create plot of each TSS group, based on mean differential expression vs. p-value.

plotted_TSS <- ggplot(results_mean_complete, aes(x = mean_diff, y = neg_log_pvalue, color = Hist_group)) +
  geom_point(size = 4) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_text(aes(x = 2.9, y = threshold, label = "0.01 Threshold"), vjust = -0.5, hjust = -1.5, color = "red") +
  theme_minimal() +
  coord_cartesian(xlim = c(2.15, 3.2)) +
  labs(title = "Comparative Analysis of Average Differential Expression \n in Target Gene Set vs Random Set", 
       x = "Average Differential Log2 Expression", 
       y = "Adjusted P-Value (-log10)") +
  theme(legend.position = "bottom",
      legend.text = element_text(size= 9),
      axis.title = element_text(size = 12),
      legend.title =element_text(hjust = 0.5,
                                 size = 10),
      plot.title = element_text(hjust = 0.5,
                                size = 14),
      panel.grid.major = element_line(color = "black",
                                      linewidth = 0.5,
                                      linetype = 1),
      panel.grid.minor = element_line(color = "black",
                                      linewidth = 0.2,
                                      linetype = 1),
      panel.background = element_rect(fill = "#FFFFF0"),
      panel.border = element_rect(fill = "transparent",
                                  color = 1,
                                  linewidth = 2)) +
  scale_color_discrete(name = "Histology Group")




# Part II. Intra-Gene-set & inter-cancer type expression analysis


# Reorder the data based on the histology group for readability
# and interpretation of graph and legend.


data_complete <- data %>% arrange(Hist_group) 
unique_Abbrvs <- unique(data_complete$TSS_Abbrvs)
data_complete$TSS_Abbrvs <- factor(data_complete$TSS_Abbrvs, levels = unique_Abbrvs)

# Q1. Boxplot showing the various distributions of RNA for each cancer type.

Bxplt_mean_RNA_by_TSS <- ggplot(data_complete, aes(x = TSS_Abbrvs, y = Mean_set, fill = Hist_group)) + 
  geom_boxplot() +
  stat_boxplot(geom = "errorbar") +
  theme_minimal() +
  labs(title = "Average Gene-set Expression by Source Tissue Type",
       y = "Log2 Expression",
       fill = "Histology Group") +
  theme(legend.position = "right",
        legend.text = element_text(size= 7),
        axis.text.x = element_text(angle= 90, 
                                   size= 8),
        axis.title.x = element_blank(),
        legend.title =element_text(hjust = 0.5,
                                   size = 11),
        plot.title = element_text(hjust = 0.5,
                                  size = 14),
        panel.grid.major = element_line(color = "black",
                                        linewidth = 0.5,
                                        linetype = 1),
        panel.grid.minor = element_line(color = "black",
                                        linewidth = 0.2,
                                        linetype = 1),
        panel.background = element_rect(fill = "#FFFFF0"),
        panel.border = element_rect(fill = "transparent",
                                    color = 1,
                                    linewidth = 2)) 



# Q1. For the ANOVA test between the groups means, first need to test for:
# 1) Normality of data within each group

LAST <- ggplot(data_complete, aes(x = Mean_set, fill = TSS_Abbrvs)) + 
  geom_histogram(alpha = 0.3, position = "identity") +
  theme_minimal() +
  labs(title = "Histogram of Gene-set Expression by Source Tissue Type",
       y = "Log2 Expression",
       fill = "TSS") +
  theme(legend.position = "right",
        legend.text = element_text(size= 7),
        axis.text.x = element_text(size= 8),
        legend.title =element_text(hjust = 0.5,
                                   size = 11),
        plot.title = element_text(hjust = 0.5,
                                  size = 14),
        panel.grid.major = element_line(color = "black",
                                        linewidth = 0.5,
                                        linetype = 1),
        panel.grid.minor = element_line(color = "black",
                                        linewidth = 0.2,
                                        linetype = 1),
        panel.background = element_rect(fill = "#FFFFF0"),
        panel.border = element_rect(fill = "transparent",
                                    color = 1,
                                    linewidth = 2)) 


ggplot(data_complete, aes(sample= Mean_set, color = TSS)) +
  stat_qq() + 
  stat_qq_line() +
  facet_wrap(~TSS) + 
  Theme(
    legend.position = "right",
    legend.text = element_text(size= 7),
    axis.text.x = element_text(size= 8),
    legend.title =element_text(hjust = 0.5,
                               size = 11),
    plot.title = element_text(hjust = 0.5,
                              size = 14),
    panel.grid.major = element_line(color = "black",
                                    linewidth = 0.5,
                                    linetype = 1),
    panel.grid.minor = element_line(color = "black",
                                    linewidth = 0.2,
                                    linetype = 1),
    panel.background = element_rect(fill = "#FFFFF0"),
    panel.border = element_rect(fill = "transparent",
                                color = 1,
                                linewidth = 2)
  ) 


# 2) Homoscedasticity of data within each group 




