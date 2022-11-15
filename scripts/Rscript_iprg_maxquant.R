##### iprg-MaxQuant

library(MSstats)
?MSstats

# First, get protein ID information
proteinGroups <- read.table("data/data_DDA_iPRG_MaxQuant/input/proteinGroups.txt", sep = "\t", header = TRUE)

# Read in MaxQuant file: evidence.txt
evi <- read.table("data/data_DDA_iPRG_MaxQuant/input/evidence.txt", sep="\t", header=TRUE)
colnames(evi)
unique(evi$Raw.file)

# Read in annotation including condition and biological replicates: annotation.csv
annot.maxquant <- read.csv("data/data_DDA_iPRG_MaxQuant/input/ABRF2015_MaxQuant_annotation.csv", header = TRUE)
annot.maxquant

?MaxQtoMSstatsFormat

# reformating and pre-processing for MaxQuant output.
# no protein with 1 peptide
input.maxquant <- MaxQtoMSstatsFormat(evidence=evi, 
                                      annotation=annot.maxquant,
                                      proteinGroups=proteinGroups,
                                      removeProtein_with1Peptide=TRUE)
head(input.maxquant)

#### Preliminary check

length(unique(input.maxquant$ProteinName)) 
sum(is.na(input.maxquant$Intensity)) 
sum(!is.na(input.maxquant$Intensity) & input.maxquant$Intensity==0)

# save your work
save(input.maxquant, file='data/data_DDA_iPRG_MaxQuant/output/input.maxquant.rda')


#load(file='output/input.maxquant.rda')

quant.maxquant <- dataProcess(raw = input.maxquant, 
                              logTrans=2, 
                              normalization = 'equalizeMedians',
                              summaryMethod = 'TMP', 
                              MBimpute=TRUE,
                              censoredInt='NA',
                             # cutoffCensored='minFeature',
                              maxQuantileforCensored = 0.999)

dataProcessPlots(data = quant.maxquant, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 'allonly',
                 address='data/data_DDA_iPRG_MaxQuant/output/ABRF_maxquant_equalMed_')

dataProcessPlots(data = quant.maxquant, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'P55249',
                 address="data/data_DDA_iPRG_MaxQuant/output/ABRF_maxquant_equalMed_P55249_")

# save your work
save(quant.maxquant, file='data/data_DDA_iPRG_MaxQuant/output/quant.maxquant.rda')


## inference
#load(file='output/quant.maxquant.rda')

unique(quant.maxquant$ProcessedData$GROUP_ORIGINAL)

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

comparison


test.maxquant <- groupComparison(contrast.matrix=comparison, data=quant.maxquant)
MaxQuant.result <- test.maxquant$ComparisonResult

# save your work
save(MaxQuant.result, file='data/data_DDA_iPRG_MaxQuant/output/MaxQuant.result.rda')
write.csv(MaxQuant.result, file='data/data_DDA_iPRG_MaxQuant/output/testResult_ABRF_maxquant.csv')

groupComparisonPlots(data = MaxQuant.result, 
                     type = 'VolcanoPlot',
                     address = 'data/data_DDA_iPRG_MaxQuant/output/testResult_ABRF_maxquant_')

groupComparisonPlots(data = MaxQuant.result, 
                     type = 'ComparisonPlot',
                     which.Protein = 'P55249',
                     address = 'data/data_DDA_iPRG_MaxQuant/output/testResult_ABRF_maxquant_')
