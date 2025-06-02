library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

selected_max_depth <- 1
selected_min_child <- 1
selected_lr <- 0.03
selected_gamma <- 0
nrounds <- 1000


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
) %>%
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
  count(freq) %>%
  as.data.frame() %>%
  spread(key = freq, value = n) %>%
  replace(is.na(.), 0)

arm_weights <- freq %>%
  mutate(total = rowSums(select(., -arm))) %>%
  mutate_at(vars(-arm, -total), list(~ 1 - round(. / total, 2))) %>%
  mutate(total = rowSums(select(., -arm, -total))) %>%
  mutate_at(vars(-arm, -total), list(~ round(. / total, 2))) %>%
  select(-total) %>%
  t() %>%
  as.data.frame() %>%
  setNames(make.unique(unlist(.[1, ]))) %>% # Convert first row to column names
  .[-1, ]

target_weights <- arm_weights[, index]
target <- as.character(colnames(arm_weights)[index])
cat("\n", target, "weights: ")
print(target_weights)
cat("\n\n")

# Manipulate the data
# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) target_weights[as.numeric(x)])
}

# Combine the data
full_df <- log_scld_tpm %>%
  merge(chr_cnv, by = 0) %>%
  column_to_rownames("Row.names") %>%
  mutate_all(as.numeric)

# Get the sample IDs from the first data frame
all_sample_ids <- rownames(full_df)

# Calculate 25% of the total samples
num_samples_to_select <- round(0.25 * length(all_sample_ids))

# Get a random sample of sample IDs
set.seed(97)
random_sample_ids <- sample(
  all_sample_ids, num_samples_to_select,
  replace = FALSE
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

xgb.train <- xgb.DMatrix(data = train_x, label = train_y, weight = weights)
xgb.test <- xgb.DMatrix(data = test_x, label = test_y)

# Train the model
# Set the parameters
xgb_params <- list(
  objective = "multi:softmax",
  eval_metric = "mlogloss",
  num_class = 3,
  max_depth = selected_max_depth,
  eta = selected_lr,
  gamma = selected_gamma,
  min_child_weight = selected_min_child
)
# Create a df to store the metrics
class_metrics <- data.frame()

# Loop over two variations
for (i in 1:2) {
  type <- if (i == 1) "Weighted" else "Unweighted"
  print(paste0("Training ", type, " model"))
  if (i == 1) {
    xgb.train <- xgb.DMatrix(data = train_x, label = train_y, weight = weights)
  } else {
    xgb.train <- xgb.DMatrix(data = train_x, label = train_y)
  }
  xgb.test <- xgb.DMatrix(data = test_x, label = test_y)

  # Train the model
  xgb_model <- xgb.train(
    params = xgb_params,
    data = xgb.train,
    watchlist = list(train = xgb.train, test = xgb.test),
    print_every_n = 10,
    nrounds = nrounds,
    early_stopping_rounds = 10
  )

  xgb_pred <- predict(xgb_model, xgb.test)
  xgb_pred <- as.data.frame(xgb_pred)
  rownames(xgb_pred) <- rownames(test_set)
  xgb_pred$actual <- test_y

  confuse_matrix <- confusionMatrix(
    factor(xgb_pred$actual, levels = c(0, 1, 2)),
    factor(xgb_pred$xgb_pred, levels = c(0, 1, 2))
  )

  # classification confusion matrix data
  new_class_metrics <- as.data.frame(confuse_matrix$byClass) %>%
    rownames_to_column("Class") %>%
    mutate(
      Class = as.integer(str_replace_all(Class, "Class: ", "")),
      Feature = target
    ) %>%
    mutate(Class = case_when(
      Class == 0 ~ "Loss",
      Class == 1 ~ "Normal",
      Class == 2 ~ "Gain"
    )) %>%
    select(-"Sensitivity", -"Pos Pred Value")

  # Gather the data into long format
  metrics_long <- gather(
    new_class_metrics,
    key = Metric,
    value = Value,
    -Class,
    -Feature
  ) %>%
    mutate(Type = type)

  # append the data
  class_metrics <- rbind(class_metrics, metrics_long)

  # Create bar plot comparing the metrics between the 3 classes
  metrics_plot <- ggplot(
    metrics_long,
    aes(x = Metric, y = Value, fill = Class)
  ) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#ABDDDE", "#CCEDB1", "#41B7C4")) +
    coord_flip() +
    labs(
      title = paste0(target, " Class Metrics"),
      subtitle = type,
      x = "Metric",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  # Save the plot
  ggsave(
    paste0(
      "Plots/",
      target,
      "_",
      type,
      "_Class_Metrics.png"
    ),
    metrics_plot,
    width = 6,
    height = 8
  )

  # Save the data
  write.csv(
    metrics_long,
    paste0("Data/Model_output/categorical/Confusion_output/",
      type, "_",
      target, "_Class_Metrics.csv"
    ),
    row.names = FALSE
  )
}
# save the data
write.csv(class_metrics,
          paste0(
            "Data/Model_output/categorical/Confusion_output/Full_Class_Metrics_",
            target,
            ".csv"
          ),
          row.names = FALSE
)

# Create bar plot comparing the metrics between the 3 classes
class_metrics_plot <- ggplot(class_metrics, aes(x = Metric, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#ABDDDE", "#CCEDB1")) +
  coord_flip() +
  labs(
    title = paste0(target, " Class Metrics"),
    x = "Metric"
  ) +
  theme_minimal() +
  theme( 
    axis.text.x = element_text(hjust = 1), 
    plot.title = element_text(hjust = 0.5)
  ) +
  facet_wrap(~Class, scales = "free_y", ncol = 3)

# Save the plot
ggsave(
  paste0(
    "Plots/",
    target,
    "Full_Class_Metrics.png"
  ),
  class_metrics_plot,
  width = 8,
  height = 6
)
