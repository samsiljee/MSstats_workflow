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

##############################
## Read Proteome Discoverer report
##############################
raw <- read.csv("data/Proteome_Discoverer_example/ControlMixture_DDA_ProteomeDiscoverer_input.csv", stringsAsFactors=F) # the data file

annot <- read.csv('data/Proteome_Discoverer_example/ControlMixture_DDA_ProteomeDiscoverer_annotation.csv')


##############################
## Make MSstats required format
##############################
quant <- PDtoMSstatsFormat(raw, 
                           annotation=annot,
                           removeProtein_with1Peptide=TRUE)

head(quant)

## count the number of proteins
length(unique(quant$ProteinName)) # 1292


##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
##############################

processed.quant <- dataProcess(quant,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               cutoffCensored="minFeature",
                               censoredInt="NA",
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

save(processed.quant, file='data/Proteome_Discoverer_example/processed.quant.rda')

##############################
## Data visualization
##############################

dataProcessPlots(processed.quant, type="QCplot", 
                 ylimDown=0, 
                 which.Protein = 'allonly',
                 width=7, height=7,  
                 address="ControlMixture_DDA_ProteomeDiscoverer_")

dataProcessPlots(processed.quant, type="Profileplot", 
                 ylimDown=0, 
                 originalPlot = TRUE,
                 summaryPlot = TRUE,
                 width=7, height=7,  
                 address="ControlMixture_DDA_ProteomeDiscoverer_")

dataProcessPlots(processed.quant, type="Conditionplot", 
                 ylimDown=0, 
                 width=7, height=7,  
                 address="ControlMixture_DDA_ProteomeDiscoverer_")


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
row.names(comparison)<-c("M1-M2", "M1-M3", "M1-M4", "M1-M5", "M2-M3", 
                         "M2-M4", "M2-M5", "M3-M4", "M3-M5", "M4-M5")


test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

##############################
## save the result
##############################

save(test.MSstats, file='test.MSstats.rda')
write.csv(test.MSstats, file='ControlMixture_DDA_ProteomeDiscoverer_testResult_byMSstats.csv')


##############################
## Visualization of result
##############################
groupComparisonPlots(data=test.MSstats, type="VolcanoPlot",
                     width=6, height=6,
                     address="ControlMixture_DDA_ProteomeDiscoverer_")

groupComparisonPlots(data=test.MSstats, type="ComparisonPlot",
                     width=6, height=6,
                     address="ControlMixture_DDA_ProteomeDiscoverer_")


