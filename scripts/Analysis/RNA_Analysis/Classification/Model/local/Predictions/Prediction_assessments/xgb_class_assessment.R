
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)


knitr::opts_knit$set(root.dir = "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

depth <- 5
min_child <- 1
lr <- 0.03
nrounds <- 1000
target <- "3p"


# Loading the test variables and data
rna_data_path <- "Data/RNA_Data/Model_Input/Train/train_"

log_scld_tpm <- read.csv(
  paste0(
    rna_data_path,
    "log_scld_tpm_soi.csv"
  ),
  row.names = 1
)

# Arm Level Aneuploidies
# Load the data
chr_cnv <- read_tsv(
  "Data/CIN_Features/CNV_Data/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
)%>% 
  replace(is.na(.), 0) %>%
  select(-"Aneuploidy Score", -"Type") %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>% 
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>% 
  mutate_all(~ replace(., . == 0, 1)) %>% 
  mutate_all(~ replace(., . == -1, 0))

# Calc the class weights
freq <- chr_cnv %>% 
  as.data.frame() %>% 
  gather(key = "arm", value = "freq") %>% 
  group_by(arm) %>%
  count(freq)%>%
  as.data.frame() %>%
  spread(key = freq, value = n) %>%
  replace(is.na(.), 0)

arm_weights <- freq %>%
  mutate(total = rowSums(select(., -arm))) %>%  
  mutate_at(vars(-arm, -total), list(~1- round(. / total,2))) %>% 
  mutate(total = rowSums(select(., -arm,-total))) %>%
  mutate_at(vars(-arm, -total), list(~ round(./total, 2))) %>%
  select(-total) %>% 
  t() %>% 
  as.data.frame() %>%
  setNames(make.unique(unlist(.[1, ]))) %>%   # Convert first row to column names
  .[-1, ]  

target_weights <- arm_weights[[target]]
# target_weights <- c(0.5, 0.5, 0.5)


# Manipulate the data
# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) target_weights[as.numeric(x)])}  

# Combine the data
full_df <- log_scld_tpm %>% 
  merge(chr_cnv, by = 0) %>% 
  column_to_rownames( "Row.names") %>%
  mutate_all(as.numeric) 

# Get the sample IDs from the first data frame
all_sample_ids <- rownames(full_df)

# Calculate 25% of the total samples
num_samples_to_select <- round(0.25 * length(all_sample_ids))

# Get a random sample of sample IDs
set.seed(97)
random_sample_ids <- sample(
  all_sample_ids, num_samples_to_select, replace = FALSE
)

# Define the train and test sets
train_set <- full_df[!rownames(full_df) %in% random_sample_ids, ]
test_set <- full_df[rownames(full_df) %in% random_sample_ids, ]

train_y_factor <- factor(train_set[, target], levels = c(0, 1, 2))
weights <- as.numeric(feature_digit_function(train_y_factor))

train_x <- as.matrix(train_set[, c(colnames(log_scld_tpm))])
train_y <- as.numeric(train_set[, target])

test_x <- as.matrix(test_set[, c(colnames(log_scld_tpm))])
test_y <- as.numeric(test_set[, target])

xgb.train = xgb.DMatrix(data=train_x,label=train_y, weight = weights)
xgb.test = xgb.DMatrix(data=test_x,label=test_y)


# Train the model
# Set the parameters 
xgb_params <- list(
  objective = "multi:softmax",
  eval_metric = "mlogloss",
  num_class = 3,
  max_depth = depth,
  eta = lr,
  gamma = 0,
  min_child_weight = min_child
)

# Train the model
xgb_model <- xgb.train(
  params = xgb_params,
  data = xgb.train,
  watchlist = list(train = xgb.train, test = xgb.test),
  print_every_n = 10,
  nrounds = nrounds,
  early_stopping_rounds = 100
)

xgb_pred <- predict(xgb_model, xgb.test)
xgb_pred <- as.data.frame(xgb_pred)
rownames(xgb_pred) <- rownames(test_set)
xgb_pred$actual <- test_y

confuse_matrix <- confusionMatrix(
  factor(xgb_pred$actual, levels = c(0, 1, 2)),
  factor(xgb_pred$xgb_pred, levels = c(0, 1, 2)))

matrix_df <- as.data.frame(as.table(confuse_matrix$table))

wide_df <- pivot_wider(matrix_df, names_from = Reference, values_from = Freq, values_fill = list(frequency = 0)) %>% 
  column_to_rownames("Prediction") %>% 
  sweep(., 2, colSums(.), "/") %>% 
  round(., 2) %>% 
  rownames_to_column("Prediction")

long_df <- pivot_longer(wide_df, cols = -Prediction, names_to = "Reference", values_to = "NormalizedFreq") %>% 
  mutate(NormalizedFreq = round(NormalizedFreq, 2))

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(long_df, aes(x = Reference, y = Prediction, fill = NormalizedFreq)) +
  geom_tile(color = "black",
            linetype = 1,
            lwd = 0.4) +  # Use white lines to separate the tiles
  geom_text(aes(label = NormalizedFreq), color = "black") +  # Add text labels
  scale_fill_gradient(low = "white", high = "red3") +  # Colors can be adjusted
  labs(title = "Normalised Multi-class Confusion Matrix",
       subtitle = paste0("Feaure: ",target),
       x = "Actual Class",
       y = "Predicted Class") +
  theme_minimal() +
  theme(axis.text.x = element_text( hjust = 1),  # Adjust text angle for better legibility
        plot.title = element_text(hjust = 0.5))  # Center the plot title



