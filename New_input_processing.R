psms_datatable <- as.data.table(PSMs)

# Full peptide-protein graphs ----
pp_orig <- unique(psms_datatable[, .(PeptideSequence = `Annotated.Sequence`, ProteinName = `Protein.Accessions`)])
pp <- cSplit(pp_orig, sep = ";", direction = "long", drop = FALSE, splitCols = "ProteinName")
pp[, ProteinName := stringr::str_replace_all(ProteinName, "\\-1", "")]
pp_graph <- createPeptideProteinGraph(pp)

# Cluster identification ----
pp <- addClusterMembership(pp, pp_graph)
input_all_prots <- merge(psms_datatable, pp, by.x = "Annotated.Sequence", by.y = "PeptideSequence", allow.cartesian = T, all.x = T, all.y = T)
input_all_prots[, `Master.Protein.Accessions` := NULL]
input_all_prots[, `Protein.Accessions` := NULL]

# MSstatsTMT pre-processing ----
rm(pp_graph, pp_orig)
gc()

procd <- MSstatsTMT::PDtoMSstatsTMTFormat(
  input_all_prots, annotations,
  which.proteinid = "ProteinName",
  useNumProteinsColumn = FALSE,
  useUniquePeptide = FALSE,
  rmPSM_withfewMea_withinRun = TRUE,
  rmProtein_with1Feature = TRUE,
  summaryforMultipleRows = max,
  use_log_file = FALSE,
  append = FALSE,
  verbose = TRUE
)

saveRDS(procd, "processed_data/tpp/procd.RDS")
procd <- readRDS("processed_data/tpp/procd.RDS")
procd <- as.data.table(procd)
rm(psms_datatable, pp, input_all_prots)
gc()

# Process isoforms ---
procd[, log2Intensity := log(Intensity, 2)]
graph_procd <- createPeptideProteinGraph(procd)
procd <- addClusterMembership(procd, graph_procd)
procd_iso <- processIsoforms(procd, T, T, F)

saveRDS(procd_iso, "processed_data/tpp/procd_iso.RDS")
procd_iso <- readRDS("processed_data/tpp/procd_iso.RDS")

# Normalization ----
procd_iso[, Intensity := 2^log2Intensity]
procd_iso <- normalizeSharedPeptides(procd_iso)
procd_iso_graph <- createPeptideProteinGraph(procd_iso)
procd_iso <- addClusterMembership(procd_iso, procd_iso_graph)
rm(procd_iso_graph)

# Cluster statistics ----
procd_iso <- getClusterStatistics(procd_iso, TRUE)

saveRDS(procd_iso, "processed_data/tpp/sub_int_cls.RDS")

# Functions ----
normalizeProteins <- function(summarized_data) {
  n_runs <- data.table::uniqueN(summarized_data$Run, na.rm = TRUE)
  if ((n_runs > 1)) {
    group_info <- unique(summarized_data$Condition)
    if (is.element("Norm", group_info)) {
      summarized_data[!is.na(Abundance), `:=`(NumRuns, data.table::uniqueN(Run,
        na.rm = TRUE
      )), by = "Protein"]
      summarized_data[!is.na(Abundance), `:=`(NumRunsWithNorm, MSstatsTMT:::.countRunsWithNorm(
        Run,
        Condition
      )), by = "Protein"]
      summarized_data[!is.na(Abundance), `:=`(
        NormalizationAbundance,
        MSstatsTMT:::.getNormalizationAbundance(Abundance, Condition)
      ),
      by = c("Protein", "Run")
      ]
      summarized_data[!is.na(Abundance), `:=`(MedianNormalized, MSstatsTMT:::.getRunsMedian(.SD)),
        by = "Protein", .SDcols = c("Run", "NormalizationAbundance")
      ]
      summarized_data[!is.na(Abundance), `:=`(Diff, MedianNormalized -
        NormalizationAbundance)]
      summarized_data[!is.na(Abundance), `:=`(
        NormalizedAbundance,
        Abundance + Diff
      )]
      summarized_data[, `:=`(Abundance, ifelse(NumRuns > 1 & NumRunsWithNorm >
        1, NormalizedAbundance, Abundance))]
      summarized_data[, `:=`(Diff, NULL)]
    } else {
      NULL
    }
  }
  summarized_data[, list(
    Mixture, TechRepMixture, Run, Channel, Protein,
    Abundance, BioReplicate, Condition
  )]
}

# Input data ----
int_cls_tbl <- readRDS("processed_data/tpp/sub_int_cls.RDS")

# Protein cluster processing and descriptive statistics ----
int_cls_each_uni <- int_cls_tbl[(HasUnique)]
int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
int_cls_each_uni[, HasUnique := any(IsUnique), by = "ProteinName"]
int_cls_each_uni <- int_cls_each_uni[(HasUnique)]
int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]

split <- split(int_cls_each_uni[NumProteins > 1], int_cls_each_uni[NumProteins > 1, Cluster])

length(split)

# Summarization ----
shared_summaries_int <- lapply(
  split,
  function(x) {
    print(unique(x$Cluster))
    tryCatch(
      {
        getWeightedProteinSummary(
          x,
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
      },
      error = function(e) NULL
    )
  }
)

unique_summaries_int <- lapply(
  split,
  function(x) {
    print(unique(x$Cluster))
    if (nrow(x[(IsUnique)]) > 0) {
      getWeightedProteinSummary(
        x,
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
    } else {
      NULL
    }
  }
)
all_summaries_int <- lapply(
  split,
  function(x) {
    print(unique(x$Cluster))
    lapply(split(x, x$ProteinName), function(y) {
      y$IsUnique <- TRUE
      getWeightedProteinSummary(
        x,
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
    })
  }
)

table(sapply(shared_summaries_int, is.null))
table(sapply(unique_summaries_int, is.null))
table(sapply(all_summaries_int, is.null))

saveRDS(shared_summaries_int, "processed_data/tpp/sh_summs_new.RDS")
saveRDS(unique_summaries_int, "processed_data/tpp/un_summs_new.RDS")
saveRDS(all_summaries_int, "processed_data/tpp/al_summs_new.RDS")

shared_summaries_int <- readRDS("processed_data/tpp/sh_summs.RDS")
unique_summaries_int <- readRDS("processed_data/tpp/un_summs.RDS")
all_summaries_int <- readRDS("processed_data/tpp/al_summs.RDS")

protein_data_shared <- rbindlist(lapply(shared_summaries_int, proteinData))
protein_data_unique <- rbindlist(lapply(unique_summaries_int, proteinData))
protein_data_all <- rbindlist(lapply(all_summaries_int, function(x) rbindlist(lapply(x, proteinData))))

feat_data_shared <- rbindlist(lapply(shared_summaries_int, featureData))
feat_data_unique <- rbindlist(lapply(unique_summaries_int, featureData))
feat_data_all <- rbindlist(lapply(all_summaries_int, function(x) rbindlist(lapply(x, featureData))))

uniqueN(feat_data_shared$ProteinName)
uniqueN(feat_data_unique$ProteinName)
uniqueN(feat_data_all$ProteinName)

# Group comparison -----
print("Done! Woo!")
contrast_matrix
gc_sh <- MSstatsTMT::groupComparisonTMT(
  list(
    ProteinLevelData = protein_data_shared,
    FeatureLevelData = feat_data_shared
  ), contrast_matrix,
  moderated = FALSE, # Steph: TRUE
  adj.method = "BH",
  remove_norm_channel = TRUE,
  remove_empty_channel = TRUE,
  save_fitted_models = FALSE,
  use_log_file = FALSE,
  verbose = TRUE
)

gc_un <- MSstatsTMT::groupComparisonTMT(list(
  ProteinLevelData = protein_data_unique,
  FeatureLevelData = feat_data_unique
), contrast_matrix,
moderated = FALSE, # Steph: TRUE
adj.method = "BH",
remove_norm_channel = TRUE,
remove_empty_channel = TRUE,
save_fitted_models = FALSE,
use_log_file = FALSE,
verbose = TRUE)

gc_al <- MSstatsTMT::groupComparisonTMT(list(
  ProteinLevelData = protein_data_all,
  FeatureLevelData = feat_data_all
), contrast_matrix,
moderated = FALSE, # Steph: TRUE
adj.method = "BH",
remove_norm_channel = TRUE,
remove_empty_channel = TRUE,
save_fitted_models = FALSE,
use_log_file = FALSE,
verbose = TRUE)

gc_sh_dt <- as.data.table(gc_sh$ComparisonResult)
gc_un_dt <- as.data.table(gc_un$ComparisonResult)
gc_al_dt <- as.data.table(gc_al$ComparisonResult)

write.csv(gc_sh_dt, file = "processed_data/tpp/MSstatsTMTResults_shared.csv")
write.csv(gc_un_dt, file = "processed_data/tpp/MSstatsTMTResults_unique.csv")
write.csv(gc_al_dt, file = "processed_data/tpp/MSstatsTMTResults_all.csv")

writeLines(
  capture.output(sessionInfo()),
  paste0(
    "processed_data/tpp/sessionInfo_",
    format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
    ".txt"
  )
)
