library(dplyr)
library(ggplot2)


path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/"
setwd(path)

tss_meta <- read.csv("Data/Other/TCGA_meta/tissueSourceSite.tsv", sep = "\t")
abbrv_meta <- read.csv("Data/Other/TCGA_meta/bcrBatchCode.tsv", sep = "\t")
meta <- left_join(
  tss_meta %>%
    dplyr::select(c("TSS.Code", "Study.Name")) %>%
    distinct() %>%
    sapply(trimws) %>%
    as.data.frame(),
  abbrv_meta %>%
    dplyr::select(c("Study.Abbreviation", "Study.Name")) %>%
    distinct() %>%
    sapply(trimws) %>%
    as.data.frame(),
  by = "Study.Name"
)

example_set <- read.csv(
  "Data/RNA_Data/Model_Input/Full/log_scld_tpm_soi.csv"
) %>%
  select("X") %>%
  mutate(other = vapply(
                        strsplit(X, "\\."),
                        `[`, 2, FUN.VALUE = character(1))
  ) %>%
  mutate(abbrvs = meta$Study.Abbreviation[match(.$other, meta$TSS.Code)])

example_set <- example_set %>%
  mutate(abbrvs = ifelse(is.na(abbrvs), tss_meta$BCR[is.na(tss_meta$TSS.Code)], abbrvs))

# ordered <- example_set %>%
#   select("abbrvs") %>%
#   mutate(group = ifelse(row_number() < 8000, "Train", "Test")) %>%
#   group_by(group, abbrvs) %>%
#   summarise(Count = n(), .groups = "drop") %>%
#   group_by(group) %>%
#   mutate(Countp = round(Count / sum(Count), 3)) %>%
#   ungroup()
# 
# 
# plt <- ggplot(ordered, aes(fill = group, x = abbrvs, y = Count)) +
#   geom_bar(stat = "identity")
# ggsave("ordered_sample_dist.pdf", width = 15, height = 5)


# # Get the total number of rows in example_set
# total_rows <- nrow(example_set)
# 
# # Calculate 70% of the total rows
# num_rows_to_select <- round(0.25 * total_rows)
# 
# # Get a random sample of row indices
# set.seed(99)
# random_indices <- sample(1:total_rows, num_rows_to_select)
# 
# combined_table <- example_set %>%
#   select("abbrvs") %>%
#   mutate(group = ifelse(!row_number() %in% random_indices, "Train", "Test")) %>%
#   group_by(group, abbrvs) %>%
#   summarise(Count = n(), .groups = "drop") %>%
#   group_by(group) %>%
#   mutate(Countp = round(Count / sum(Count), 3)) %>%
#   ungroup()
# 
# 
# plt2 <- ggplot(combined_table, aes(fill = group, x = abbrvs, y = Count)) +
#   geom_bar(stat = "identity")
# ggsave("sample_dist.pdf", width = 15, height = 5)

