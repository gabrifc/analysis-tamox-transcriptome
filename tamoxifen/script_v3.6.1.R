##----------------------------------------------------------------------------
## Libraries, files
##----------------------------------------------------------------------------
library("tximport")
library('DESeq2')
library('dplyr')
library('tidyr')
library('tibble')
library('BiocParallel')
library('ggplot2')
library('ggrepel')
library('eulerr')


##----------------------------------------------------------------------------
## Helper functions through the script 
##----------------------------------------------------------------------------
downloadMapping <- function(vectorOfIDs) {
  library("biomaRt")
  # print("Annotating using BioMaRt")
  zfishMart <- useMart("ensembl", dataset="drerio_gene_ensembl")
  txMapping <- getBM(attributes = c('ensembl_gene_id', 
                                    'description', 
                                    'external_gene_name',
                                    'entrezgene',
                                    'zfin_id_id'), 
                     filters = 'ensembl_gene_id',
                     values = vectorOfIDs, 
                     mart = zfishMart)
  # Clean blank mappings
  txMapping <- txMapping[!txMapping$ensembl_gene_id=="",]
  return(txMapping)
}
fixFDR <- function(resSkewed) {
  library('fdrtool')
  resSkewed <- resSkewed[ !is.na(resSkewed$padj), ]
  resSkewed <- resSkewed[ !is.na(resSkewed$pvalue), ]
  resSkewed <- resSkewed[ !is.na(resSkewed$stat), ]
  resSkewed <- resSkewed[, -which(names(resSkewed) == "padj")]
  
  FDR.resSkewed <- fdrtool(resSkewed$stat, statistic= "normal", plot = FALSE)
  FDR.resSkewed$param[1, "sd"]
  # Plot the corrected histogram of p-values
  par(mfrow=c(1, 2))
  hist(FDR.resSkewed$pval, breaks = 0:20/20,
       col = "grey50", border = "white", xlab = "CORRECTED p-values")
  hist(FDR.resSkewed$pval[resSkewed$baseMean>1], breaks = 0:20/20,
       col = "grey50", border = "white", xlab = "CORRECTED p-values")
  par(mfrow=c(1, 1))
  
  resSkewed[,"padj"]  <- p.adjust(FDR.resSkewed$pval, method = "BH")
  return(resSkewed)
}
annotateRes <- function(res) {
  resMapped <- merge(as.data.frame(res), txMapping,
                     by.x = 'row.names',
                     by.y = "ensembl_gene_id")
  resMapped <- deleteDuplicatesDataFrame(resMapped, 'Row.names')
  rownames(resMapped) <- resMapped[,1]
  resMapped <- resMapped[,-1]
}
downloadGeneLengthData <- function(id) {
  # ahDb <- query(ah, pattern = c("Danio Rerio", "GRCz11", "EnsDb")) to get the ID
  # Input: the ID of the EnsDB object.
  # Output: A vector with the gene Length Data.
  # Update the gene lengths
  library(ensembldb)
  library(AnnotationHub)
  ah <- AnnotationHub()
  ahDb <- query(ah, "EnsDb")
  EnsDb <- ahDb[[id]]
  return(lengthOf(EnsDb, of="gene"))
}
runGoSeq <- function(genes, geneLengthData, goseqMapping, goSeqMethod) {
  library('goseq')
  # Analysis
  # Fitting the Probability Weighting Function (PWF) to check if DE genes are
  # biased by lenght
  pwf=nullp(genes, "danRer11", "ensGene", bias.data = geneLengthData)
  
  # For every case below, we can only select the GO that wer are interested using
  # test.cats=c("GO:MF"): GO.MF=goseq(pwf,"danRer10","ensGene", test.cats=c("GO:MF"))
  
  # Analysis
  if (goSeqMethod == 'Wallenius') {
    GO = goseq(pwf, "danRer11", "ensGene", gene2cat = goMapping)
  } else if (goSeqMethod == 'RandomSampling') {
    GO = goseq(pwf, "danRer11", "ensGene", method="Sampling", repcnt=1000)
  } else if (goSeqMethod == 'NoBias') {
    GO = goseq(pwf, "danRer11", "ensGene", method="Hypergeometric")
  }
  
  # Output
  output <- list()
  GO$over_represented_padj <- p.adjust(GO$over_represented_pvalue, method="BH")
  GO$under_represented_padj <- p.adjust(GO$under_represented_pvalue, method="BH")
  goResults <- subset(GO, over_represented_pvalue<.05)
  #goResults <- subset(GO, over_represented_pvalue<.05)
  #goResults$padj <- p.adjust(GO$over_represented_pvalue, method="BH")
  rownames(goResults) <- goResults$category
  unEnrichedGo <- subset(GO, under_represented_pvalue<.05)
  #unEnrichedGo$padj <- p.adjust(GO$under_represented_pvalue, method="BH")
  #unEnrichedGo <- subset(GO, p.adjust(under_represented_pvalue, method="BH")<.05)
  rownames(unEnrichedGo) <- unEnrichedGo$category
  
  output$enrichedGo <- goResults
  output$unEnrichedGo <- unEnrichedGo
  return(output)
}
deleteDuplicatesDataFrame <- function(df, col) {
  dup.idx <- which(duplicated(df[col]))
  return(df[-dup.idx,])
}
deleteDuplicatesVector <- function(vector) {
  vector <- as.vector(vector)
  dup.idx = which(duplicated(vector))
  if (length(dup.idx) > 0) {
    return(vector[-dup.idx])
  } else {
    return(vector)
  }
}
renderPathways <- function(listOfPathways, foldChangeData, suffix, keggDataFolder) {
  # Function renderPathways.
  # Input: list of Pathways and gene expression
  # Output: none. Creates images for every pathway in the list.
  ## Detach dplyr because select causes error.
  detach("package:dplyr", unload=TRUE)
  library("pathview")
  # Limit for log2FC colors
  limit = list(gene=0.5, cpd=1) 
  # Low expression color
  low = list(gene="deepskyblue3",cpd="deepskyblue3")
  # Mid Expression color
  mid = list(gene="gray",cpd="gray")
  # High Expression color
  high = list(gene="gold1",cpd="gold1")
  pv.out.list <- sapply(listOfPathways, function(pid) pathview(
    gene.data = foldChangeData, pathway.id = pid, species = "dre", 
    limit=limit, low=low, mid=mid, high=high, min.nnodes = 0, node.sum="mean", 
    kegg.native = TRUE, res = 600, same.layer = FALSE,
    bins = list(gene = 20, cpd = 10), cex = 0.15, 
    out.suffix = suffix, kegg.dir = keggDataFolder))
  library('dplyr')
}
makeVolcanoPlot <- function(df = df, 
                            log2FCColumn = 'log2FoldChange',
                            sigColumn = 'svalue',
                            sigCutoff = 0.005,
                            sigLabel = 'sValue') {
  
  df <- as.data.frame(df)
  names(df)[names(df) == log2FCColumn] <- 'log2FoldChange'
  names(df)[names(df) == sigColumn] <- 'dfSigValues'
  
  dfPlotData <- as.data.frame(
    dplyr::mutate(as.data.frame(df), 
                  sig = ifelse(df$dfSigValues < sigCutoff, "Sig", "Not Sig")), 
    row.names=rownames(df))
  
  vPlot <- ggplot2::ggplot(dfPlotData, 
                           ggplot2::aes(log2FoldChange, 
                                        -log10(dfSigValues))) +
    ggplot2::geom_point(ggplot2::aes(col = sig)) +
    ggplot2::scale_color_manual(values = c("black", "red")) +
    ggplot2::ggtitle("Volcano Plot of DESeq2 analysis") +
    ggplot2::ylab(paste0('-log10(',sigLabel,')'))
  
  return(vPlot)
}

##----------------------------------------------------------------------------
## Load Data 
##----------------------------------------------------------------------------

## Mapping (we can create these files with downloadMapping)
tx2gene <- read.csv("data/tx2gene.csv", 
                    row.names=1, sep=";", stringsAsFactors=FALSE)

txMapping <- read.csv("data/txMapping.csv", 
                      row.names=1, sep=";", stringsAsFactors=FALSE)

## Update GO terms from Ensembl Biomart.
goTermsFile <- 'data/zfishGO.txt'
goMapping <- read.delim(goTermsFile, na.strings = "", stringsAsFactors = FALSE)
colnames(goMapping) <- c('Ensembl', 'GO')
goMapping <- goMapping[complete.cases(goMapping),]

# Get zfish 2 human orthology via zfin
zfish2human <- read.csv("data/zfin.txt", skip = 1,
                        sep="\t", stringsAsFactors=FALSE) %>%
  deleteDuplicatesDataFrame(df = ., 
                            col = "ZFIN.ID") %>%
  dplyr::select(ZFIN.ID, Human.Symbol, Gene.ID)

##----------------------------------------------------------------------------
## Import Salmon data with tximport 
##----------------------------------------------------------------------------
## Create the tx2gene needed for importing the runs.

# ## Load all TranscriptIDs
# all <- read.delim("C:/Users/gabri/Desktop/ralfRNAseq/all.txt")
# 
# universeIDs <- deleteDuplicatesVector(all$Name)
# 
# remove(all)
# 
# ## Create the tx2gene matrix
# tx2gene <- downloadMapping(universeIDs)
# 
# ## Save it for later runs
# write.csv2(tx2gene, file = 'tx2gene.csv')

## Load the counts
countDirectory <- "Counts/"
samples <- list.files(path = countDirectory)
## Clean the names
names(samples) <- sub(".txt", "", samples)
## Add the directory to the count names
samples <- paste0(countDirectory, samples)
samples
txi <- tximport(samples, type = "salmon", tx2gene = tx2gene)


##----------------------------------------------------------------------------
## QC 
##----------------------------------------------------------------------------

## We are not removing any sample, but the 'batch' or 'individual' effect seems
## to be greater than the treatment or infection effect (view plots), so we will
## have to take into account when designing the matrix with deseq2.

## Final PCA Plots
# rld <- rlog(ddsTxi)
# plotPCA(rld, intgroup = c('treatment', 'infection'))
# plotPCA(rld, intgroup = c('genotype'))


##----------------------------------------------------------------------------
## Load DESeq2 
##----------------------------------------------------------------------------
sampleGrouping <- read.csv('data/sampleGrouping.csv', sep=";", row.names = 1)

## Reorder sample grouping
sampleGrouping <- sampleGrouping[order(rownames(sampleGrouping)),]
all(rownames(sampleGrouping) == colnames(txi$counts))

##----------------------------------------------------------------------------
## -- Model DESeq2: Ungrouped Modelling
##----------------------------------------------------------------------------

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sampleGrouping,
                                   design = ~ genotype + treatment +infection + treatment:infection)

## Specify the controls.
ddsTxi$treatment <- relevel(ddsTxi$treatment, ref = "DMSO")
ddsTxi$infection <- relevel(ddsTxi$infection, ref = "No")

##----------------------------------------------------------------------------
## -- Check fitType, local or parametric? -> Local 
##----------------------------------------------------------------------------

## See https://support.bioconductor.org/p/81094/

# ddsPara <- DESeq(ddsTxi, fitType = "parametric")
# ddsLocal <- DESeq(ddsTxi, fitType = "local")
# 
# par(mfrow=c(1, 2))
# plotDispEsts(ddsPara, main = 'Parametric Fit')
# plotDispEsts(ddsLocal, main = 'Local Fit')
# par(mfrow=c(1, 1))
# 
# ## Calculate the ratios of median absolute residue
# logParaDispGeneEst=log(mcols(ddsPara)$dispGeneEst)
# logParaDispFit=log(mcols(ddsPara)$dispFit)
# paraRatio <- abs(median(na.omit(logParaDispGeneEst))) -
#   abs(median(na.omit(logParaDispFit)))
# 
# logLocalDispGeneEst=log(mcols(ddsLocal)$dispGeneEst)
# logLocalDispFit=log(mcols(ddsLocal)$dispFit)
# localRatio <- abs(median(na.omit(logLocalDispGeneEst))) -
#   abs(median(na.omit(logLocalDispFit)))
# 
# localRatio < paraRatio

## As the local fit "fits" better than the parametric, I'll choose local.

##----------------------------------------------------------------------------
## -- Run Deseq2
##----------------------------------------------------------------------------
dds <- DESeq(ddsTxi, fitType = "local")
resultsNames(dds)

##----------------------------------------------------------------------------
## TAM vs DMSO - No infected (results similar to Grouped modelling)
##----------------------------------------------------------------------------
resTam <- lfcShrink(dds,
                    coef = 6,
                    type = 'apeglm',
                    svalue=TRUE,
                    parallel=TRUE,
                    BPPARAM=SnowParam(4))

resTamFDR <- results(dds, 
                     name = 'treatment_TAM_vs_DMSO')
summary(resTam)
plot(resTamFDR$padj, resTam$svalue, col="blue",
     xlab="DESeq2 padj", ylab="apeglm svalue",
     xlim=c(0,1), ylim=c(0,0.5))

plot(resTamFDR$padj, resTam$svalue, col="blue",
     xlab="DESeq2 padj", ylab="apeglm svalue",
     xlim=c(0,0.15), ylim=c(0,0.005))

resTamMapped <- annotateRes(resTam)
vPlotTam <- makeVolcanoPlot(resTamMapped)
vPlotTam

write.table(resTamMapped,
           file = 'TAM_all_vs_DMSO_all_lfcSh.tsv',
           sep = "\t")

##----------------------------------------------------------------------------
## Interaction Term -> Is the infection different in Tam treated?
##----------------------------------------------------------------------------
resIntTamMm <- lfcShrink(dds,
                      coef = 9,
                      type = 'apeglm',
                      svalue=TRUE,
                      parallel=TRUE,
                      BPPARAM=SnowParam(4))

resIntTamMmFDR <- results(dds, 
                       name = 'treatmentTAM.infectionMm')
summary(resIntTamMm)
plot(resIntTamMmFDR$padj, resIntTamMm$svalue, col="blue",
     xlab="DESeq2 padj", ylab="apeglm svalue",
     xlim=c(0,1), ylim=c(0,0.5))

plot(resIntTamMmFDR$padj, resIntTamMm$svalue, col="blue",
     xlab="DESeq2 padj", ylab="apeglm svalue",
     xlim=c(0,0.25), ylim=c(0,0.005))

resIntTamMmMapped <- annotateRes(resIntTamMm)
vPlotIntTamMm <- makeVolcanoPlot(resIntTamMmMapped)
vPlotIntTamMm

# write.csv2(resIntTamMmMapped, 
#            file = 'Interaction_TAM_Mm_lfcSh.csv')

##----------------------------------------------------------------------------
## Infection Mm vs DMSO
##----------------------------------------------------------------------------
resMm <- lfcShrink(dds,
                   coef = 7,
                   type = 'apeglm',
                   svalue = TRUE)
summary(resMm)

resMmFDR <- results(dds, 
                          name = 'infection_Mm_vs_No')

resMmMapped <- annotateRes(resMm)
vPlotMm <- makeVolcanoPlot(resMmFDR,
                           sigColumn = 'padj',
                           sigCutoff = 0.05,
                           sigLabel = 'padj')
vPlotMm

# write.csv2(resMmMapped, 
#            file = 'DMSO_Mm_vs_DMSO_uninfected_lfcSh.csv')

##----------------------------------------------------------------------------
## Infection TamMm vs Tam
##----------------------------------------------------------------------------

# The infection effect for Tam
# This is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in Tam compared to DMSO).
resTamMm <- results(dds, 
                    list(c("infection_Mm_vs_No","treatmentTAM.infectionMm")))

resTamMm <- lfcShrink(dds,
                   contrast = list(c("infection_Mm_vs_No","treatmentTAM.infectionMm")),
                   type = 'apeglm',
                   svalue = TRUE)
summary(resMm)
resMmMapped <- annotateRes(resMm)
vPlotMm <- makeVolcanoPlot(resMmMapped)
vPlotMm

##----------------------------------------------------------------------------
## Gene Ontology Enrichment
##----------------------------------------------------------------------------

## Download geneLengthData Needed for GOseq from the latest Ensembl Release
#geneLengthData <- downloadGeneLengthData('AH64906')

testGenes <- resIntTamMm
sValueCutOff <- 0.005
fileNameComparison <- 'Interaction_TAM_Mm_lfcSh'

## Update the geneLengthData with our rnaseq genes
universeIDs <- rownames(testGenes)
genesWithLengthData <- universeIDs[universeIDs %in% names(geneLengthData)]
geneLengthData <- geneLengthData[names(geneLengthData) %in% genesWithLengthData]

## Change for each comparison: 
## Get the significant genes
significantGenes <- as.data.frame(testGenes) %>%
  rownames_to_column('gene') %>%
  dplyr::filter(svalue < sValueCutOff) %>%
  column_to_rownames('gene')

genes <- as.integer(genesWithLengthData %in% rownames(significantGenes))
names(genes) <- genesWithLengthData
table(genes)
## runGoSeq
goSeqOutput <- runGoSeq(genes = genes, 
                        geneLengthData = geneLengthData, 
                        goseqMapping =  goMapping, 
                        goSeqMethod = 'Wallenius')
goResults <- goSeqOutput$enrichedGo
unEnrichedGO <- goSeqOutput$unEnrichedGo

## write the enrichment results
# write.csv2(as.data.frame(goResults), 
#            file=paste0(fileNameComparison, ".goseq.GO.csv"))

# write.csv2(as.data.frame(unEnrichedGO),
#            file=paste0(fileNameComparison, ".unEnrichedGO.csv"))

##----------------------------------------------------------------------------
## Pathway Analysis
##----------------------------------------------------------------------------

comparison <- resTam
sValueCutOff <- 0.005
fileName <- 'TAM_uninfected_vs_DMSO_uninfected'

## For the pathway analysis we will need EntrezIDs. We can get that from biomaRt
resMapped <- annotateRes(comparison)

resMappedSig <- as.data.frame(resMapped) %>%
  rownames_to_column('gene') %>%
  dplyr::filter(svalue < sValueCutOff) %>%
  column_to_rownames('gene')

# Get Entrez ID of significant genes and all the dataset
significantGenes <- na.omit(resMappedSig$entrezgene)
significantGenes <- deleteDuplicatesVector(significantGenes)

universeGenes <- na.omit(resMapped$entrezgene)
universeGenes <- deleteDuplicatesVector(universeGenes)

library('limma')

pathwayResults <- kegga(significantGenes, universe = universeGenes, 
                        species = "Dr")

pathwayResults$padj <- p.adjust(pathwayResults$P.DE, method="BH")

#pathwayResults <- subset(pathwayResults, P.DE<.1)

pathwayResults <- subset(pathwayResults, p.adjust(P.DE, method="BH")<.1)

pathwayFileName <- paste0(fileName, ".kegga.csv")
write.csv2(pathwayResults, file = pathwayFileName)

# Visualization variables 
foldChangeData <- resMapped[,c("entrezgene", "log2FoldChange")]
foldChangeData <- deleteDuplicatesDataFrame(foldChangeData, "entrezgene")
foldChangeList <- foldChangeData[,2]
names(foldChangeList) <- foldChangeData$entrezgene

# pathwayResults, extract the number of the pathway. e.g. from 
# "path:dre00022" to 00022
listOfPathways <- c(substr(c(rownames(pathwayResults)), 9, 13))
suffix <- fileName
keggDataFolder <- "pathwayRenders/keggData"

# renderPathways(listOfPathways, foldChangeData = foldChangeList, 
#                suffix = suffix, keggDataFolder = keggDataFolder)

##----------------------------------------------------------------------------
## GSEA - Create Ranks
##----------------------------------------------------------------------------

comparison <- resTam
fileName <- 'TAM_uninfected_vs_DMSO_uninfected'

resMapped <- annotateRes(comparison)
df <- data.frame(resMapped[,c("zfin_id_id", "log2FoldChange", "svalue")])
df <- df[order(df$svalue),]
df <- deleteDuplicatesDataFrame(df = df, col = "zfin_id_id")
#df <- deleteDuplicatesDataFrame(df = df, col = "Homolog")

#df$rank <- -log10(df$padj)*df$log2FC
df$rank <- -log10(df$svalue)*sign(df$log2FoldChange)

rankedDF <- df[,c("zfin_id_id", "rank")]
rankedDF <- rankedDF[complete.cases(rankedDF),]
#rankedDF$Symbol <- toupper(rankedDF$Symbol)

outputFileName <- paste0(fileName, ".geneSet.svalueSign.zfin.rnk")

write.table(rankedDF, file = outputFileName, sep = "\t", 
            row.names=FALSE, col.names = FALSE, quote = FALSE)

##----------------------------------------------------------------------------
## Some plots :)
##----------------------------------------------------------------------------

resIntTamMmSig <- resIntTamMmMapped %>%
  rownames_to_column('gene') %>%
  dplyr::filter(svalue < 0.005) %>%
  select(gene, external_gene_name)

resIntTamMmSigGenes <- resIntTamMmSig %>%
  select(gene) %>%
  unlist()

remove(countData)

countData <- lapply(resIntTamMmSigGenes, function(x) plotCounts(dds, x, c("treatment", "infection"),returnData = TRUE))
for(i in 1:length(countData)) countData[[i]]$gene <- resIntTamMmSig[i, 'external_gene_name']

countData <- do.call(rbind, countData) %>%
  dplyr::filter(treatment != 'AMIO') %>%
  dplyr::mutate(infectionNames = ifelse(infection == 'No', "PBS", "Mm"))
  

countData$infectionNames <- factor(countData$infectionNames, 
                                   levels=c("PBS", "Mm"))

ggplot(countData, 
       aes(x = infectionNames, 
           y = count, 
           colour = treatment)) + 
  #scale_y_log10() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), 
             size = 2) + 
  stat_summary(aes(y = count,
                   group=treatment,
                   colour = treatment), 
               fun.y=mean,
               geom="line",
               size = 1) +
  facet_wrap(~gene, ncol = 7, scales = "free_y") +
  ggtitle(label =  'Genes whose expression during Mycobacterium marinum infection depends on the treatment',
          subtitle = 'Significance: svalue < 0.005; Dots represent normalized read counts per library; Lines connect the means in each treatment.') +
  xlab('') +
  ylab('Gene Read Counts') + 
  theme_bw() +
  theme(legend.position="bottom") 

plotCounts(dds, gene = "ENSDARG00000005972", c("treatment", "infection"))


##----------------------------------------------------------------------------
## EulerR
##----------------------------------------------------------------------------
resTamSig <- resTamMapped %>%
  rownames_to_column('gene') %>%
  dplyr::filter(svalue < 0.005) %>%
  select(gene) %>%
  unlist() %>%
  unname()
resMmSig <- resMmMapped %>%
  rownames_to_column('gene') %>%
  dplyr::filter(svalue < 0.005) %>%
  select(gene) %>%
  unlist() %>%
  unname()

eulerList <- list(Tamoxifen = resTamSig, 
                  Infection = resMmSig)

plot(euler(eulerList))
euler(eulerList)

common <- resMmSig[resMmSig %in% resTamSig]
commonNames <- resTamMapped[common, "external_gene_name"]
commonNames

commonGenes <- cbind(resTamMapped[common, ], resMmMapped[common, ])





commonCountData <- lapply(common, function(x) plotCounts(dds, x, c("treatment", "infection"),returnData = TRUE))
for(i in 1:length(commonCountData)) commonCountData[[i]]$gene <- commonGenes[i, 'external_gene_name']

commonCountData <- do.call(rbind, commonCountData) %>%
  dplyr::filter(treatment != 'AMIO') %>%
  dplyr::mutate(treatmentNames = paste(treatment, infection, sep = '_'))


commonCountData$treatmentNames <- factor(commonCountData$treatmentNames, 
                                   levels=c("DMSO_No", "DMSO_Mm", 
                                            "TAM_No", "TAM_Mm"))

ggplot(commonCountData, 
       aes(x = treatmentNames, 
           y = count, 
           colour = treatment)) + 
  #scale_y_log10() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), 
             size = 2) + 
  stat_summary(aes(y = count,
                   group=treatment,
                   colour = treatment), 
               fun.y=mean,
               geom="line",
               size = 1) +
  facet_wrap(~gene, ncol = 7, scales = "free_y") +
  ggtitle(label =  'Genes significantly regulated by both tamoxifen and Mm infection',
          subtitle = 'Significance: svalue < 0.005; Dots represent normalized read counts per library; Lines connect the means in each treatment.') +
  xlab('') +
  ylab('Gene Read Counts') + 
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45,hjust=1)) 
