# Load test datasets
load("~/Documents/Coding/MSstats_workflow/test_input/Input_crash.rda")
load("~/Documents/Coding/MSstats_workflow/test_input/Input_error.rda")
load("~/Documents/Coding/MSstats_workflow/test_input/Input_fine.rda")

# Code to check uniqueness of rows in data
rows_fine <- nrow(Input_fine)
rows_error <- nrow(Input_error)
rows_crash <- nrow(Input_crash)

columns_to_select <- c(
  "Cluster",
  "ProteinName",
  "PeptideSequence",
  "Run",
  "Channel",
  "Charge",
  "PSM",
  "Intensity",
  "log2Intensity",
  "log2IntensityNormalized",
  "TechRepMixture",
  "BioReplicate",
  "Condition",
  "NumProteins",
  "NumPeptides",
  "TotalSize",
  "NumProteinsPerPeptide",
  "NumPeptidesPerProtein",
  "IsUnique",
  "HasUnique",
  "AnyHasUnique",
  "EachHasUnique"
)

columns_to_select <- c(
  "TechRepMixture",
  "BioReplicate",
  "Condition",
  "Intensity",
  "log2Intensity",
  "log2IntensityNormalized",
  "ProteinName",
  "PSM"
)

summarise_duplicates <- function(Input_data) {
  Input_summary <- Input_data %>% mutate(Index = paste(
    TechRepMixture,
    BioReplicate,
    Condition,
    Intensity,
    log2Intensity,
    log2IntensityNormalized,
    ProteinName,
    PSM, sep = "_")) %>%
    group_by(Index) %>%
    summarise(Count = n()) %>%
    filter(Count > 1) %>%
    arrange(desc(Count))
  
  return(Input_summary)
}

# Finding duplicate rows
Input_fine_summary <- summarise_duplicates(Input_fine)
Input_error_summary <- summarise_duplicates(Input_error)
Input_crash_summary <- summarise_duplicates(Input_crash)

# When filtering by  ProteinName, PSM, TechRepMixture, BioReplicate, Condition, Intensity, log2Intensity, and log2IntensityNormalized, the duplicated rows all come in duplicates.
# It's worth noting that the duplicated rows all have "NA" values for all intesity values (and derived values).
# There are no such duplicates in the erroring data.
# It's also worth noting that the duplicates only occur in TMT channels 126 and 135N, both channels that I used for the pooled reference.

rows_fine - nrow(unique(select(Input_fine, columns_to_select)))
rows_error - nrow(unique(select(Input_error, columns_to_select)))
rows_crash - nrow(unique(select(Input_crash, columns_to_select)))

# Try removing channels causing errors in crashing data
# Doesn't cause R to crash, but still stops with the same error as the error dataset.
MSstats_summarised <- getWeightedProteinSummary(
  filter(Input_crash, !Channel %in% c("126", "135N")),
  norm = "p_norm",
  norm_parameter = 1,
  weights_mode = "contributions",
  tolerance = 0.1,
  max_iter = 10,
  initial_summary = "unique",
  weights_penalty = FALSE,
  weights_penalty_param = 0.1,
  save_weights_history = FALSE,
  save_convergence_history = FALSE
)

# Filter out proteins causing the crash
Crash_proteins <- filter(Input_crash_index, Index %in% Input_crash_summary$Index) %>% .$ProteinName

# Try again without these proteins
MSstats_summarised <- getWeightedProteinSummary(
  filter(Input_crash, !ProteinName %in% Crash_proteins),
  norm = "p_norm",
  norm_parameter = 1,
  weights_mode = "contributions",
  tolerance = 0.1,
  max_iter = 10,
  initial_summary = "unique",
  weights_penalty = FALSE,
  weights_penalty_param = 0.1,
  save_weights_history = FALSE,
  save_convergence_history = FALSE
)

# # Get cluster sizes
# cluster_sizes <- MSstatsTMT_input %>%
#   group_by(Cluster) %>%
#   summarise(Cluster_count = n())
# 
# # Initiate blank index
# cluster_iteration_index <- numeric()
# 
# # Loop through and select a cluster to use as an example for each one
# for(i in sort(unique(cluster_sizes$Cluster_count))) {
#   cluster_iteration_index <- c(
#     cluster_iteration_index,
#     as.numeric(dplyr::filter(cluster_sizes, Cluster_count == i)[1,1])
#   )
# }
# 
# # Manually change some clusters not functioning
# cluster_iteration_index[cluster_iteration_index == 420] <- 495
# cluster_iteration_index[cluster_iteration_index == 772] <- 1182
# 
# # Remove clusters already tested
# cluster_iteration_index <- cluster_iteration_index[48:length(cluster_iteration_index)]
# 
# # testing loop
# index_count <- 0

# Make a list of cluster exclusions
stop_exclusion_list <- c(8, 12, 18, 24, 25, 30)
crash_exclusion_list <- c(1, 22, 32)

indexes <- unique(MSstatsTMT_input$Cluster)[!unique(MSstatsTMT_input$Cluster) %in% 1:max(c(stop_exclusion_list, crash_exclusion_list))]

for(i in indexes) {
  print(paste("Now testing cluster:", i))
  MSstats_summarised <- getWeightedProteinSummary(
    dplyr::filter(MSstatsTMT_input, Cluster == i),
    norm = "p_norm",
    norm_parameter = 1,
    weights_mode = "contributions",
    tolerance = 0.1,
    max_iter = 10,
    initial_summary = "unique",
    weights_penalty = FALSE,
    weights_penalty_param = 0.1,
    save_weights_history = FALSE,
    save_convergence_history = FALSE
  )
}