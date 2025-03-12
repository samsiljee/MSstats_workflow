# Libraries ----
library(data.table)
library(MSstatsWeightedSummary)
library(ggplot2)
library(splitstackshape)
library(parallel)
library(ggVennDiagram)
library(patchwork)
# Raw data ----
# ## onePot
psms = fread("PSM_subset.csv")
load("./annotations.rda")
# Full peptide-protein graphs ----
# ## onePot
pp_orig = unique(psms[, .(PeptideSequence = `Annotated.Sequence`, ProteinName = `Protein.Accessions`)])
pp = cSplit(pp_orig, sep = "; ", direction = "long", drop = FALSE, splitCols = "ProteinName")
pp[, ProteinName := stringr::str_replace_all(ProteinName, "\\-1", "")]
pp_graph = createPeptideProteinGraph(pp)
# # Cluster identification ----
# ## onePot
pp = addClusterMembership(pp, pp_graph)
input_all_prots = merge(psms, pp, by.x = "Annotated.Sequence", by.y = "PeptideSequence", allow.cartesian = T, all.x = T, all.y = T)
input_all_prots[, `Master.Protein.Accessions` := NULL]
input_all_prots[, `Protein.Accessions` := NULL]
# # MSstatsTMT pre-processing ----
## onePot
procd = MSstatsTMT::PDtoMSstatsTMTFormat(input_all_prots, annotations,
                                         which.proteinid = "ProteinName",
                                         useNumProteinsColumn = FALSE,
                                         useUniquePeptide = FALSE,
                                         rmPSM_withfewMea_withinRun = FALSE,
                                         rmProtein_with1Feature = FALSE,
                                         summaryforMultipleRows = max,
                                         use_log_file = FALSE,
                                         append=FALSE,
                                         verbose=TRUE)
procd = as.data.table(procd)
# Process isoforms ---
# ## onePot
procd[, log2Intensity := log(Intensity, 2)]
graph_procd = createPeptideProteinGraph(procd)
procd = addClusterMembership(procd, graph_procd)
procd[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]

ggplot(procd, 
       aes(x = Channel, y = log2Intensity, group = PSM, color = IsUnique)) +
  geom_line() +
  theme_bw() +
  facet_grid(~ProteinName) +
  theme_bw()

procd_merge_single_identical = processIsoforms(procd, T, T, F)
procd_merge_single_identical[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
ggplot(procd_merge_single_identical, 
       aes(x = Channel, y = log2Intensity, group = PSM, color = IsUnique)) +
  geom_line() +
  theme_bw() +
  facet_grid(~ProteinName) +
  theme_bw()

# Optional steps: 
# merging or filtering isoforms if there is not enough feature-level information 
# procd_iso = processIsoforms(procd, T, T, F)
# Normalization (feature level, median-based)
# procd_iso[, Intensity := 2 ^ log2Intensity]
# procd_iso = normalizeSharedPeptides(procd_iso)
# procd_iso_graph = createPeptideProteinGraph(procd_iso)
# procd_iso = addClusterMembership(procd_iso, procd_iso_graph)
# Cluster statistics 
# procd_iso = getClusterStatistics(procd_iso, TRUE)
