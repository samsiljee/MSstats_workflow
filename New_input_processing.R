# New test to import the data - this time using the vignette: https://vitek-lab.github.io/MSstatsWeightedSummary/articles/shared_workflow.html

# Required columns ProteinName, PeptideSequence, Charge, PSM, Run, Channel, Intensity, Condition, BioReplicate, Mixture, TechRepMixture

library(tidyr) # Used for separate_rows function

# Format using original MSstatsTMT processing - will also deal with fractions
MSstatsTMT_input <- MSstatsTMT::PDtoMSstatsTMTFormat(
  input = PSMs[1:1000,],
  annotation = annotations,
  which.proteinid = "Protein.Accessions", # Default: Protein Accessions
  useNumProteinsColumn = FALSE,
  useUniquePeptide = FALSE,
  rmPSM_withfewMea_withinRun = FALSE, # Default: TRUE
  rmProtein_with1Feature = FALSE,
  summaryforMultipleRows = sum,
  use_log_file = TRUE,
  append = FALSE,
  verbose = TRUE,
  log_file_path = paste0(
    output_dir,
    "/logs/PDtoMSstatsTMTFormat_log_",
    format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
    ".txt"
  )
)

psms <- as.data.table(MSstatsTMT_input)

# Normalisation
psms <- normalizeSharedPeptides(psms)

# Add protein group clusters
pp_graph <- createPeptideProteinGraph(psms)
psms <- addClusterMembership(psms, pp_graph)

# Load PSMs from PD output
psms <- PSMs %>%
  filter(Spectrum.File == "Mix_TMT_F1_2_dot_5ul.raw") %>%
  # Pivot from wide to long format
  pivot_longer(
    cols = `Abundance:.126`:`Abundance:.135N`,
    names_to = "Channel",
    names_prefix = "Abundance:.",
    values_to = "Intensity"
  ) %>%
  # Rename columns and add column to merge with annotations
  mutate(
    ProteinName = Protein.Accessions,
    PeptideSequence = Annotated.Sequence,
    PSM = paste(Annotated.Sequence, Charge, sep = "_"),
    key = paste0(Spectrum.File, Channel)
  ) %>%
  # Merge with annotations
  merge(
    mutate(annotations, key = paste0(Run, Channel)),
    all.x = TRUE,
    by = "key"
  ) %>%
  mutate(Channel = Channel.y) %>%
  # Select required columns
  select(ProteinName, PeptideSequence, Charge, PSM, Run, Channel, Intensity, Condition, BioReplicate, Mixture, TechRepMixture) %>%
  separate_rows(ProteinName, sep = "; ") %>%
  # Convert to data.table for MSstatsWeightedSummary
  as.data.table()

# Add protein group clusters
pp_graph <- createPeptideProteinGraph(psms)
psms <- addClusterMembership(psms, pp_graph)

# Process isoforms
psms <- processIsoforms(psms, T, T, F)

# Normalisation
psms <- normalizeSharedPeptides(psms)

# Add protein group clusters
pp_graph <- createPeptideProteinGraph(psms)
psms <- addClusterMembership(psms, pp_graph)

# Get summaries
MSstats_summarised <- getWeightedProteinSummary(
  psms,
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
