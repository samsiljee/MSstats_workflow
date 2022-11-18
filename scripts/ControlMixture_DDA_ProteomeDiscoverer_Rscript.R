# This script was not in the original tutorial of the video, however it was referred to at the end
# The data can be found at; https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=5dfbaf705b894d369aaed6f60d51000e#%7B%22table_sort_history%22%3A%22main.collection_asc%22%7D
# The details of the original study can be found at; https://massive.ucsd.edu/ProteoSAFe/reanalysis.jsp?task=5dfbaf705b894d369aaed6f60d51000e
# I have changed some of the path names to reflect the new paths

##############################
#### Analysis in MSstats
##############################

##############################
## Load MSstats package
##############################
library(MSstats)
library(tidyverse)
##############################
## Read Proteome Discoverer report
##############################
raw <- read.csv("data/Proteome_Discoverer_example/input/ControlMixture_DDA_ProteomeDiscoverer_input.csv", stringsAsFactors = TRUE)

annot <- read.csv('data/Proteome_Discoverer_example/input/ControlMixture_DDA_ProteomeDiscoverer_annotation.csv', stringsAsFactors = TRUE)

##############################
## Make MSstats required format
##############################
input <- PDtoMSstatsFormat(raw, 
                           annotation=annot,
                           removeProtein_with1Peptide=FALSE)

head(input)

## count the number of proteins
input$ProteinName %>% unique %>% length # 1292


##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
## Argument deleted, to be compatible with newer version of `MSstats`; 'cutoffCensored = "minFeature"'
## See discussion; https://groups.google.com/g/msstats/c/WLZQz1-QdiA/m/3yjqu4TBAAAJ
##############################

processed.quant <- dataProcess(input,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               censoredInt="NA",
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

save(processed.quant, file='data/Proteome_Discoverer_example/output/processed.quant.rda')

##############################
## Data visualization
##############################

dataProcessPlots(processed.quant, type="QCplot", 
                 ylimDown=0, 
                 which.Protein = 'allonly',
                 width=7, height=7,  
                 address="data/Proteome_Discoverer_example/output/ControlMixture_DDA_ProteomeDiscoverer_")

dataProcessPlots(processed.quant, type="Profileplot", 
                 ylimDown=0, 
                 originalPlot = TRUE,
                 summaryPlot = TRUE,
                 width=7, height=7,  
                 address="data/Proteome_Discoverer_example/output/ControlMixture_DDA_ProteomeDiscoverer_")

dataProcessPlots(processed.quant, type="Conditionplot", 
                 ylimDown=0, 
                 width=7, height=7,  
                 address="data/Proteome_Discoverer_example/output/ControlMixture_DDA_ProteomeDiscoverer_")


##############################
## Model-based comparison + adjust p-value
##############################

comparison1<-matrix(c(1,-1,0,0,0),nrow=1)
comparison2<-matrix(c(1,0,-1,0,0),nrow=1)
comparison3<-matrix(c(1,0,0,-1,0),nrow=1)
comparison4<-matrix(c(1,0,0,0,-1),nrow=1)
comparison5<-matrix(c(0,1,-1,0,0),nrow=1)
comparison6<-matrix(c(0,1,0,-1,0),nrow=1)
comparison7<-matrix(c(0,1,0,0,-1),nrow=1)
comparison8<-matrix(c(0,0,1,-1,0),nrow=1)
comparison9<-matrix(c(0,0,1,0,-1),nrow=1)
comparison10<-matrix(c(0,0,0,1,-1),nrow=1)
comparison<-rbind(comparison1,comparison2, comparison3, comparison4, comparison5, 
                  comparison6, comparison7, comparison8, comparison9, comparison10)
row.names(comparison)<-c("Condition1-Condition2", "Condition1-Condition3", "Condition1-Condition4", "Condition1-Condition5", "Condition2-Condition3", 
                         "Condition2-Condition4", "Condition2-Condition5", "Condition3-Condition4", "Condition3-Condition5", "Condition4-Condition5")
colnames(comparison) <- paste("Condition", 1:5, sep = "")

test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

##############################
## save the result
##############################

save(test.MSstats, file='data/Proteome_Discoverer_example/output/test.MSstats.rda')
write.csv(test.MSstats, file='data/Proteome_Discoverer_example/output/ControlMixture_DDA_ProteomeDiscoverer_testResult_byMSstats.csv')


##############################
## Visualization of result
##############################
groupComparisonPlots(data=test.MSstats, type="VolcanoPlot",
                     width=6, height=6,
                     address="data/Proteome_Discoverer_example/output/ControlMixture_DDA_ProteomeDiscoverer_")

groupComparisonPlots(data=test.MSstats, type="ComparisonPlot",
                     width=6, height=6,
                     address="data/Proteome_Discoverer_example/output/ControlMixture_DDA_ProteomeDiscoverer_")


