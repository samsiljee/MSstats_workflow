##### iprg-Skyline

library(MSstats)
?MSstats

# Read output from skyline 
raw.skyline <- read.csv(file="data/data_DDA_iPRG_Skyline/input/ABRF2015_Skyline_report.csv")

# Check the first 6 rows of dataset
head(raw.skyline)
# total number of unique protein name
length(unique(raw.skyline$Protein))

# several isotopic peaks for peptide charge
unique(raw.skyline$FragmentIon)

# unique FileName, which is MS run.
unique(raw.skyline$FileName)

# 'Truncated' column
unique(raw.skyline$Truncated)

# count table for 'Truncated' column
xtabs(~Truncated, raw.skyline)

# count which 'Truncated' is 'True'
sum(raw.skyline$Truncated == 'True')


## read annotation file
annot.skyline <- read.csv(file="input/ABRF2015_Skyline_annotation.csv")
annot.skyline

setdiff(unique(raw.skyline$FileName), annot.skyline$Run)
setdiff(annot.skyline$Run, unique(raw.skyline$FileName))

## bad example
annot.wrong <- read.csv(file="input/ABRF2015_Skyline_annotation_wrong_example.csv")
annot.wrong

setdiff(unique(raw.skyline$FileName), annot.wrong$Run)
setdiff(annot.wrong$Run, unique(raw.skyline$FileName))


# reformating and pre-processing for Skyline output.
input.skyline <- SkylinetoMSstatsFormat(raw.skyline, 
                                        annotation=annot.skyline,
                                        removeProtein_with1Feature = TRUE)
head(input.skyline)
unique(input.skyline$Run)

## Preliminary check

length(unique(input.skyline$ProteinName)) 
sum(is.na(input.skyline$Intensity)) 
sum(!is.na(input.skyline$Intensity) & input.skyline$Intensity==0)

# save the work
save(input.skyline, file='output/input.skyline.rda')

## Load the pre-processed data of Skyline output
load(file='output/input.skyline.rda')

## data processing
quant.skyline <- dataProcess(raw = input.skyline, 
                             logTrans=2, 
                             normalization = 'equalizeMedians',
                             summaryMethod = 'TMP', 
                             MBimpute=TRUE,
                             censoredInt='0', ## important for Skyline
                             cutoffCensored='minFeature',
                             maxQuantileforCensored = 0.999)

save(quant.skyline, file='output/quant.skyline.rda')

# show the name of outputs
names(quant.skyline)

# show reformated and normalized data.
# 'ABUNDANCE' column has normalized log2 transformed intensities.
head(quant.skyline$ProcessedData)

# This table includes run-level summarized log2 intensities. (column : LogIntensities)
# Now one summarized log2 intensities per Protein and Run.
# NumMeasuredFeature : show how many features are used for run-level summarization.
#         If there is no missing value, it should be the number of features in certain protein.
# MissingPercentage : the number of missing features / the number of features in certain protein.
head(quant.skyline$RunlevelData)

# show which summarization method is used.
head(quant.skyline$SummaryMethod)

# QC plot for normalized data with equalize median method
dataProcessPlots(data = quant.skyline, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 'allonly',
                 address='output/ABRF_skyline_equalizeNorm_')

# if you have many MS runs, adjust width of plot (makd wider)
# Profile plot for the data with equalized median method
dataProcessPlots(data = quant.skyline, 
                 type="Profileplot", 
                 width=7, height=7,
                 address="output/ABRF_skyline_equalizeNorm_")

dataProcessPlots(data = quant.skyline, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P44015|VAC2_YEAST',
                 address="output/ABRF_skyline_equalizeNorm_P44015_")

dataProcessPlots(data = quant.skyline, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P55249|ZRT4_YEAST',
                 address="output/ABRF_skyline_equalizeNorm_P55249_")

# not run
dataProcessPlots(data = quant.skyline, 
                 type="conditionplot", 
                 width=7, height=7,
                 address="output/ABRF_skyline_equalizeNorm_")

dataProcessPlots(data = quant.skyline, 
                 type="conditionplot", 
                 width=7, height=7,
                 which.Protein = 'sp|P44015|VAC2_YEAST',
                 address="output/ABRF_skyline_equalizeNorm_P44015_")

#### No imputation, TMP summarization only
quant.skyline.TMPonly <- dataProcess(raw = input.skyline, 
                                     logTrans=2, 
                                     summaryMethod = 'TMP', 
                                     MBimpute=FALSE, ##
                                     censoredInt='0',
                                     cutoffCensored='minFeature',
                                     maxQuantileforCensored = 0.999)

#### feature selection
quant.skyline.ftrslct <- dataProcess(raw = input.skyline, 
                                     logTrans=2, 
                                     summaryMethod = 'TMP', 
                                     MBimpute=TRUE,
                                     censoredInt='0',
                                     cutoffCensored='minFeature',
                                     maxQuantileforCensored = 0.999,
                                     featureSubset="highQuality",
                                     remove_uninformative_feature_outlier=TRUE)

## inference
#load(file='output/quant.skyline.rda')

unique(quant.skyline$ProcessedData$GROUP_ORIGINAL)

comparison1 <- matrix(c(-1,1,0,0),nrow=1)
comparison2 <- matrix(c(-1,0,1,0),nrow=1)
comparison3 <- matrix(c(-1,0,0,1),nrow=1)
comparison4 <- matrix(c(0,-1,1,0),nrow=1)
comparison5 <- matrix(c(0,-1,0,1),nrow=1)
comparison6 <- matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison) <- c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

comparison

test.skyline <- groupComparison(contrast.matrix=comparison, data=quant.skyline)


#### Save the comparison result 

# Let's save the testing result as rdata and .csv file.

Skyline.result <- test.skyline$ComparisonResult

save(Skyline.result, file='output/Skyline.result.rda')
write.csv(Skyline.result, file='output/testResult_ABRF_skyline.csv')


#### subset of significant comparisons
#Let's inspect the results to see what proteins are changing significantly between Diseased and Healthy.

head(Skyline.result)

# select subset of rows with adj.pvalue < 0.05
SignificantProteins <- 
    Skyline.result[Skyline.result$adj.pvalue < 0.05, ]

nrow(SignificantProteins)

# select subset of rows with adj.pvalue < 0.05 and log2FC > 2
SignificantProteinsUpInDiseased <- SignificantProteins[SignificantProteins$log2FC > 2 ,]

nrow(SignificantProteinsUpInDiseased)

### Visualization of differentially abundant proteins
groupComparisonPlots(data = Skyline.result, 
                     type = 'VolcanoPlot',
                     address = 'output/testResult_ABRF_skyline_')

groupComparisonPlots(data = Skyline.result, 
                     type = 'VolcanoPlot',
                     sig = 0.05, 
                     FCcutoff = 2^2, 
                     address = 'output/testResult_ABRF_skyline_FCcutoff4_')

groupComparisonPlots(Skyline.result, 
                     type="ComparisonPlot", 
                     address="output/testResult_ABRF_skyline_")

Skyline.result[Skyline.result$Protein == 'sp|P44015|VAC2_YEAST', ]


#### Calculating statistical power 
?designSampleSize

# calculate the power
test.power <- designSampleSize(data = test.skyline$fittedmodel, 
                               desiredFC = c(1.1, 1.6), 
                               FDR = 0.05,
                               power = TRUE,
                               numSample = 3)
test.power

#### Visualizing the relationship between desired fold-change and power
designSampleSizePlots(data = test.power)

#### Designing sample size for desired fold-change
# Minimal number of biological replicates per condition
samplesize <- designSampleSize(data = test.skyline$fittedmodel, 
                               desiredFC = c(1.1, 1.6), 
                               FDR = 0.05,
                               power = 0.9,
                               numSample = TRUE)
samplesize


#### Visualizing the relationship between desired fold-change and mininum sample size number
designSampleSizePlots(data = samplesize)


### Protein subject quantification 
?quantification

sampleQuant <- quantification(quant.skyline)
head(sampleQuant)

groupQuant <- quantification(quant.skyline, type='Group')
head(groupQuant)
