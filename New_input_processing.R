psms_onepot <- as.data.table(PSMs[1:1000, ])

# Full peptide-protein graphs ----
## onePot
onepot_pp_orig <- unique(psms_onepot[, .(PeptideSequence = `Annotated.Sequence`, ProteinName = `Protein.Accessions`)])
onepot_pp <- cSplit(onepot_pp_orig, sep = ";", direction = "long", drop = FALSE, splitCols = "ProteinName")
onepot_pp[, ProteinName := stringr::str_replace_all(ProteinName, "\\-1", "")]
onepot_pp_graph <- createPeptideProteinGraph(onepot_pp)

# Cluster identification ----
## onePot
onepot_pp <- addClusterMembership(onepot_pp, onepot_pp_graph)
onepot_input_all_prots <- merge(psms_onepot, onepot_pp, by.x = "Annotated.Sequence", by.y = "PeptideSequence", allow.cartesian = T, all.x = T, all.y = T)
onepot_input_all_prots[, `Master.Protein.Accessions` := NULL]
onepot_input_all_prots[, `Protein.Accessions` := NULL]

# MSstatsTMT pre-processing ----
rm(onepot_pp_graph, curve_pp_graph, onepot_pp_orig, curve_pp_raw)
gc()
# onePot
onepot_procd <- MSstatsTMT::PDtoMSstatsTMTFormat(
  onepot_input_all_prots, annotations,
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

saveRDS(onepot_procd, "processed_data/onepot_tpp/onepot_procd.RDS")
onepot_procd <- readRDS("processed_data/onepot_tpp/onepot_procd.RDS")
onepot_procd <- as.data.table(onepot_procd)
rm(psms_onepot, onepot_pp, onepot_input_all_prots, curve_data_pp, curve_pp, curve_input)
gc()

# Process isoforms ---
## onePot
onepot_procd[, log2Intensity := log(Intensity, 2)]
onepot_graph_procd <- createPeptideProteinGraph(onepot_procd)
onepot_procd <- addClusterMembership(onepot_procd, onepot_graph_procd)
onepot_procd_iso <- processIsoforms(onepot_procd, T, T, F)

saveRDS(onepot_procd_iso, "processed_data/onepot_tpp/onepot_procd_iso.RDS")
onepot_procd_iso <- readRDS("processed_data/onepot_tpp/onepot_procd_iso.RDS")

# Normalization ----
onepot_procd_iso[, Intensity := 2^log2Intensity]
onepot_procd_iso <- normalizeSharedPeptides(onepot_procd_iso)
onepot_procd_iso_graph <- createPeptideProteinGraph(onepot_procd_iso)
onepot_procd_iso <- addClusterMembership(onepot_procd_iso, onepot_procd_iso_graph)
rm(onepot_procd_iso_graph)

# Cluster statistics ----
onepot_procd_iso <- getClusterStatistics(onepot_procd_iso, TRUE)

saveRDS(onepot_procd_iso, "processed_data/onepot_tpp/onepot_sub_int_cls.RDS")

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
onepot_int_cls_tbl <- readRDS("processed_data/onepot_tpp/onepot_sub_int_cls.RDS")

# Protein cluster processing and descriptive statistics ----
onepot_int_cls_each_uni <- onepot_int_cls_tbl[(HasUnique)]
onepot_int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
onepot_int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]
onepot_int_cls_each_uni[, HasUnique := any(IsUnique), by = "ProteinName"]
onepot_int_cls_each_uni <- onepot_int_cls_each_uni[(HasUnique)]
onepot_int_cls_each_uni[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
onepot_int_cls_each_uni[, NumProteins := uniqueN(ProteinName), by = "Cluster"]

onepot_split <- split(onepot_int_cls_each_uni[NumProteins > 1], onepot_int_cls_each_uni[NumProteins > 1, Cluster])

length(onepot_split)

# Summarization ----
onepot_shared_summaries_int <- lapply(
  onepot_split,
  function(x) {
    print(unique(x$Cluster))
    tryCatch(
      {
        getWeightedProteinSummary(x, "Huber", 1e-6,
          max_iter = 100, tolerance = 1e-2, initial_summary = "flat"
        )
      },
      error = function(e) NULL
    )
  }
)

onepot_unique_summaries_int <- lapply(
  onepot_split,
  function(x) {
    print(unique(x$Cluster))
    if (nrow(x[(IsUnique)]) > 0) {
      getWeightedProteinSummary(x[(IsUnique)], "Huber", 1e-6, max_iter = 100, tolerance = 1e-2)
    } else {
      NULL
    }
  }
)
onepot_all_summaries_int <- lapply(
  onepot_split,
  function(x) {
    print(unique(x$Cluster))
    lapply(split(x, x$ProteinName), function(y) {
      y$IsUnique <- TRUE
      getWeightedProteinSummary(y, "Huber", 1e-6, max_iter = 100, tolerance = 1e-2)
    })
  }
)

table(sapply(onepot_shared_summaries_int, is.null))
table(sapply(onepot_unique_summaries_int, is.null))
table(sapply(onepot_all_summaries_int, is.null))

saveRDS(onepot_shared_summaries_int, "processed_data/onepot_tpp/onepot_sh_summs.RDS")
saveRDS(onepot_unique_summaries_int, "processed_data/onepot_tpp/onepot_un_summs.RDS")
saveRDS(onepot_all_summaries_int, "processed_data/onepot_tpp/onepot_al_summs.RDS")

onepot_shared_summaries_int <- readRDS("processed_data/onepot_tpp/onepot_sh_summs.RDS")
onepot_unique_summaries_int <- readRDS("processed_data/onepot_tpp/onepot_un_summs.RDS")
onepot_all_summaries_int <- readRDS("processed_data/onepot_tpp/onepot_al_summs.RDS")

onepot_protein_data_shared <- rbindlist(lapply(onepot_shared_summaries_int, proteinData))
onepot_protein_data_unique <- rbindlist(lapply(onepot_unique_summaries_int, proteinData))
onepot_protein_data_all <- rbindlist(lapply(onepot_all_summaries_int, function(x) rbindlist(lapply(x, proteinData))))

onepot_feat_data_shared <- rbindlist(lapply(onepot_shared_summaries_int, featureData))
onepot_feat_data_unique <- rbindlist(lapply(onepot_unique_summaries_int, featureData))
onepot_feat_data_all <- rbindlist(lapply(onepot_all_summaries_int, function(x) rbindlist(lapply(x, featureData))))

uniqueN(onepot_feat_data_shared$ProteinName)
uniqueN(onepot_feat_data_unique$ProteinName)
uniqueN(onepot_feat_data_all$ProteinName)

# Group comparison -----
cm_onepot <- contrast_matrix
gc_sh_onepot <- MSstatsTMT::groupComparisonTMT(
  list(
    ProteinLevelData = onepot_protein_data_shared,
    FeatureLevelData = onepot_feat_data_shared
  ), cm_onepot,
  use_log_file = FALSE
)
gc_un_onepot <- MSstatsTMT::groupComparisonTMT(list(
  ProteinLevelData = onepot_protein_data_unique,
  FeatureLevelData = onepot_feat_data_unique
), cm_onepot, use_log_file = FALSE)
gc_al_onepot <- MSstatsTMT::groupComparisonTMT(list(
  ProteinLevelData = onepot_protein_data_all,
  FeatureLevelData = onepot_feat_data_all
), cm_onepot, use_log_file = FALSE)

gc_sh_dt_onepot <- as.data.table(gc_sh_onepot$ComparisonResult)
gc_un_dt_onepot <- as.data.table(gc_un_onepot$ComparisonResult)
gc_al_dt_onepot <- as.data.table(gc_al_onepot$ComparisonResult)
