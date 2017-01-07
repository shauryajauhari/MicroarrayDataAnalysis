# Installing GEOquery

source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")

# Loading GEO file with GEOquery

library(Biobase)
library(GEOquery)

#Download GPL file, put it in the current directory, and load it:
gpl570 <- getGEO('GPL570', destdir=".")

#Or, open an existing GPL file:
gpl570 <- getGEO(filename='GPL570.soft')

# Handpicked description(three columns: ID, Gene Symbol, Gene Title).

Table(gpl570) [c("ID","Gene Symbol","Gene Title")]
IDs <- attr(dataTable(gpl570), "table")[, c("ID", "Gene Symbol", "Gene Title")]


# Extract the expression values from the dataset
# line 64 contains field names
DS_Main <- read.table("GSE26910_series_matrix.txt.gz", skip = 63, header = TRUE, sep = "\t", fill = TRUE)

# Remove the last line from the matrix that says "!series_matrix_table_end"
DS_Main <- DS_Main[-54676, ]

# Merging the annotation information to the expression values matrix and rejecting null entries.

names(IDs)[1] <- "ID_REF"
DS <- merge(IDs,DS_Main, by = "ID_REF")
DS[DS == ""] <- NA
DS <- na.omit(DS)

# Reordering of respective breast cancer and prostate cancer datsets.
# Prostate Normal [1:6], Prostate Tumor [7:12], Breast Normal [13:18], Breast Tumor [19:24]

WorkDS <- DS [c(4,6,8,10,12,14, 5,7,9,11,13,15, 16,18,20,22,24,26, 17,19,21,23,25,27)]

# RowMeans calculation

ProstateNormalMean <- rowMeans(log2(WorkDS[,1:6]))
ProstateTumorMean <- rowMeans(log2(WorkDS[,7:12]))
BreastNormalMean <- rowMeans(log2(WorkDS[,13:18]))
BreastTumorMean <- rowMeans(log2(WorkDS[,19:24]))

# MA-Plot

par(mfrow=c(1,2))
ProstateMean <- rowMeans(log2(WorkDS[, 1:12]))
BreastMean <- rowMeans(log2(WorkDS[, 13:24]))
plot(ProstateMean, ProstateTumorMean-ProstateNormalMean, main="MA Plot for Prostate data", pch=16, cex=0.35)
hold()
plot(BreastMean, BreastTumorMean-BreastNormalMean, main="MA Plot for Breast data", pch=16, cex=0.35)

# Rough draft of extreme probes

DS[which.min(BreastTumorMean-BreastNormalMean), ] ### most negatively expressed breast gene
DS[which.min(ProstateTumorMean-ProstateNormalMean), ] ### most negatively expressed prostate gene
DS[which.max(ProstateTumorMean-ProstateNormalMean), ] ### most positively expressed prostate gene
DS[which.max(BreastTumorMean-BreastNormalMean), ] ### most positively expressed breast gene


# Standard Deviation calculation for t-test

install.packages(genefilter)
library(genefilter)
ProstateNormalSD <- rowSds(log2(WorkDS[,1:6]))
ProstateTumorSD <- rowSds(log2(WorkDS[,7:12]))
BreastNormalSD <- rowSds(log2(WorkDS[,13:18]))
BreastTumorSD <- rowSds(log2(WorkDS[,19:24]))

# t-test calculation and histogram plot

par(mfrow=c(1,2))
Prostate_ttest <- (ProstateTumorMean-ProstateNormalMean)/sqrt(ProstateTumorSD^2/6 + ProstateNormalSD^2/6)
hist(Prostate_ttest,nclass=100)
hold()
Breast_ttest <- (BreastTumorMean-BreastNormalMean)/sqrt(BreastTumorSD^2/6 + BreastNormalSD^2/6)
hist(Breast_ttest, nclass=100)

# p-value calculation and histogram plot

Prostate_pval <- 2*(1-pt(abs(Prostate_ttest),5))
Breast_pval <- 2*(1-pt(abs(Breast_ttest),5))
par(mfrow=c(1,2))
hist(Prostate_pval, nclass=100)
hold()
hist(Breast_pval, nclass = 100)

# volcano Plot

par(mfrow=c(1,2))
plot(ProstateTumorMean-ProstateNormalMean, -log10(Prostate_pval), main ="Volcano Plot@Prostate tissue", xlab= "Sample Mean Difference", ylab= "-log10(p value)", pch=16, cex=0.35)
hold()
plot(BreastTumorMean-BreastNormalMean, -log10(Breast_pval), main ="Volcano Plot@Breast tissue", xlab= "Sample Mean Difference", ylab= "-log10(p value)", pch=16, cex=0.35)


# Boxplots for the normal data and its log transformed version.(Log2 transformation applied)

par(mfrow = c(1, 2))
boxplot(WorkDS, col = c(2,3,2,3,2,3,2,3,2,3,2,3), main = "Expression values pre-normalization", 
        xlab = "Slides", ylab = "Expression", las = 2, cex.axis = 0.7)
hold()
boxplot(log2(WorkDS), col = c(2,3,2,3,2,3,2,3,2,3,2,3), main = "Expression values post-log-transformation", 
        xlab = "Slides", ylab = "Expression", las = 2, cex.axis = 0.7)
abline(0, 0, col = "black")


# Check Normality

par(mfrow=c(1,2))
qqnorm(Prostate_ttest, main = "QQ Plot@Prostate Data")
qqline(Prostate_ttest)

hold()

qqnorm(Breast_ttest, main = "QQ Plot@Breast Data")
qqline(Breast_ttest)


# Elucidating genes with particular p-values.

for (i in c(0.01, 0.05, 0.001, 1e-04, 1e-05, 1e-06, 1e-07)) 
  print(paste("genes with p-values smaller than",i, length(which(Prostate_pval < i))))
for (i in c(0.01, 0.05, 0.001, 1e-04, 1e-05, 1e-06, 1e-07)) 
  print(paste("genes with p-values smaller than",i, length(which(Breast_pval < i))))


# Plot heatmap of differentially expressed genes: Genes are differentially expressed if its p-value is under a given threshold, which must be smaller than the usual 0.05 or 0.01 due to multiplicity of tests

BreastDEGenes <- data.frame(which(Breast_pval < 0.01))
ProstateDEGenes <- data.frame(which(Prostate_pval < 0.01))
ProstateDEGenesData <- ProstateDEGenes[ ,1]
BreastDEGenesData <- BreastDEGenes[ ,1]

ProstateData <- as.matrix(WorkDS[ProstateDEGenesData, 1:12])
heatmap(ProstateData, col = topo.colors(100), cexRow = 0.5)

BreastData <- as.matrix(WorkDS[BreastDEGenesData, 13:24])
heatmap(BreastData, col = topo.colors(100), cexRow = 0.5)


# List of differentially expressed genes.

#Breast Data
BDEG <- matrix(nrow = nrow(BreastDEGenes), ncol = 1)
for(i in 1:nrow(BreastDEGenes))  BDEG[i,]<- paste(DS[BreastDEGenes[i,], "ID_REF"])
BDEG <- as.data.frame(BDEG)
names(BDEG)[1] <- "ID_REF"
FinalBDEG <- merge(BDEG,DS)
BDEG <- merge(BDEG, IDs, by = 'ID_REF')
view(BDEG)

#Prostate Data
PDEG <- matrix(nrow = nrow(ProstateDEGenes),ncol = 1)
for(i in 1:nrow(ProstateDEGenes))  PDEG[i,] <- paste(DS[ProstateDEGenes[i,], "ID_REF"])
PDEG <- as.data.frame(PDEG)
names(PDEG)[1] <- "ID_REF"
FinalPDEG <- merge(PDEG,DS)
PDEG <- merge(PDEG, IDs, by = 'ID_REF')
view(PDEG)


##Intersecting transcripts in breast and prostate cancer types as marked in the dataset

BDEG$match <- match(BDEG$location, PDEG$location, nomatch=0)


# Reordering of respective breast cancer and prostate cancer datsets.
# Prostate Normal [1:6], Prostate Tumor [7:12], Breast Normal [13:18], Breast Tumor [19:24]

FinalPDEG <- FinalPDEG [c(4,6,8,10,12,14, 5,7,9,11,13,15, 16,18,20,22,24,26, 17,19,21,23,25,27)]
WorkFinalPDEG <- FinalPDEG[1:12]


FinalBDEG <- FinalBDEG [c(4,6,8,10,12,14, 5,7,9,11,13,15, 16,18,20,22,24,26, 17,19,21,23,25,27)]
WorkFinalBDEG <- FinalBDEG[13:24]


##Prostate data regrerssion analysis(linear model)

par(mfrow=c(1,6))
plot(log2(WorkFinalPDEG$GSM662756),log2(WorkFinalPDEG$GSM662757),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662756) ~ log2(WorkFinalPDEG$GSM662757)), col= 1)
plot(log2(WorkFinalPDEG$GSM662758),log2(WorkFinalPDEG$GSM662759),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662758) ~ log2(WorkFinalPDEG$GSM662759)), col= 1)
plot(log2(WorkFinalPDEG$GSM662760),log2(WorkFinalPDEG$GSM662761),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662760) ~ log2(WorkFinalPDEG$GSM662761)), col= 1)
plot(log2(WorkFinalPDEG$GSM662762),log2(WorkFinalPDEG$GSM662763),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662762) ~ log2(WorkFinalPDEG$GSM662763)), col= 1)
plot(log2(WorkFinalPDEG$GSM662764),log2(WorkFinalPDEG$GSM662765),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662764) ~ log2(WorkFinalPDEG$GSM662765)), col= 1)
plot(log2(WorkFinalPDEG$GSM662766),log2(WorkFinalPDEG$GSM662767),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662766) ~ log2(WorkFinalPDEG$GSM662767)), col= 1)


##Gene Set Enrichment Analysis

library(genefilter)
library(GSEABase)
Breast_GSEA <- GeneSetCollection(WorkFinalBDEG, setType = KEGGCollection())
Prostate_GSEA <- GeneSetCollection(WorkFinalPDEG, setType = KEGGCollection())


##Correlation Analysis
##Breast

WorkFinalBDEG <- read.csv("WorkFinalBDEG_GSEAFiltered.csv") ## Import filtered annotation file from MeV.
btemp <- WorkFinalBDEG
btemp$ID_REF <- NULL
btemp <- log2(btemp)
pairs(btemp)
BreastCorrelationMatrix <- cor(t(as.matrix(btemp)))
BreastCorMat <- as.data.frame(BreastCorrelationMatrix)
rownames(BreastCorMat) <- WorkFinalBDEG$ID_REF
colnames(BreastCorMat) <- WorkFinalBDEG$ID_REF

##Prostate

WorkFinalPDEG <- read.csv("WorkFinalPDEG_GSEAFiltered.csv") ## Import filtered annotation file from MeV.
ptemp <- WorkFinalPDEG
ptemp$ID_REF <- NULL
ptemp <- log2(ptemp)
pairs(ptemp)
ProstateCorrelationMatrix <- cor(t(as.matrix(ptemp)))
ProCorMat<- as.data.frame(ProstateCorrelationMatrix)
rownames(ProCorMat) <- WorkFinalPDEG$ID_REF
colnames(ProCorMat)<- WorkFinalPDEG$ID_REF


### Feature Selection: Clustering of robustly entwined genes.

install.packages("gplots")
install.packages("Hmisc")
library(Hmisc)
library(gplots)

heatmap.2(ProstateCorrelationMatrix, main="Hierarchical Cluster", dendrogram="column",trace="none",col=greenred(10))
heatmap.2(1-abs(ProstateCorrelationMatrix), distfun=as.dist, trace="none")
heatmap.2(BreastCorrelationMatrix, main="Hierarchical Cluster", dendrogram="column",trace="none",col=greenred(10))
heatmap.2(1-abs(BreastCorrelationMatrix), distfun=as.dist, trace="none")

##Prostate Data

library(caret)
HighlyCorrelated <- findCorrelation(ProstateCorrelationMatrix, cutoff = 0.95, verbose = TRUE, names = FALSE)
print(HighlyCorrelated)
WorkFinalPDEG[HighlyCorrelated,1]
IDs[WorkFinalPDEG[HighlyCorrelated,1],c(2,3)]


for(i in 2:nrow(BreastCorMat))
{
  for(j in 1:ncol(BreastCorMat)-1)
  {
     if(i>j)
      {
      out <- c (rownames(BreastCorMat[i,]), colnames(BreastCorMat[j]), BreastCorMat[i,j])
      write.table(out, file="output.txt", append=TRUE, sep= " ")
      }
     else
     break
  }
}


#Network Ready Matrix Format (Function) // Credit: http://www.sthda.com

flattenCorrMatrix <- function(cmat, pmat) {
  ut <- upper.tri(cmat)
  data.frame(
    row = rownames(cmat)[row(cmat)[ut]],
    column = rownames(cmat)[col(cmat)[ut]],
    cor  =(cmat)[ut],
    p = pmat[ut]
  )
}

library(Hmisc)

btemp <- as.matrix(btemp)
rownames(btemp)<- WorkFinalBDEG$ID_REF
BreastNet <-rcorr(t(btemp))
BreastNetworkInputMatrix<- flattenCorrMatrix(BreastNet$r, BreastNet$P)

#lets map the gene names to row and column entries
BreastNetworkInputMatrix$row <-IDs[WorkFinalBDEG[BreastNetworkInputMatrix$row,1],2]
BreastNetworkInputMatrix$column <-IDs[WorkFinalBDEG[BreastNetworkInputMatrix$column,1],2]


ptemp <- as.matrix(ptemp)
rownames(ptemp)<- WorkFinalPDEG$ID_REF
ProstateNet <-rcorr(t(ptemp))
ProstateNetworkInputMatrix <- flattenCorrMatrix(ProstateNet$r, ProstateNet$P)
ProstateNetworkInputMatrix$row <-IDs[WorkFinalPDEG[ProstateNetworkInputMatrix$row,1],2]
ProstateNetworkInputMatrix$column <-IDs[WorkFinalPDEG[ProstateNetworkInputMatrix$column,1],2]

symnum(BreastCorrelationMatrix)
symnum(ProstateCorrelationMatrix)

install.packages("corrplot")
library(corrplot)
corrplot(BreastCorrelationMatrix, type="upper", order="hclust", tl.col="black", tl.srt=45)
corrplot(ProstateCorrelationMatrix, type="upper", order="hclust", tl.col="black", tl.srt=45)

install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(BreastCorrelationMatrix, histogram= TRUE, pch= 19)
chart.Correlation(ProstateCorrelationMatrix, histogram= TRUE, pch= 19)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = BreastCorrelationMatrix, col = col, symm = TRUE)
heatmap(x = ProstateCorrelationMatrix, col = col, symm = TRUE)

## Optimize network ready correlation and p-values matrix
## top candidates which manifest low p-value and high correlation.

BreastFinal <- BreastNetworkInputMatrix[which(abs(BreastNetworkInputMatrix$cor) > 0.95 | BreastNetworkInputMatrix$p < 0.0000001), c(1,2,3,4)]
ProstateFinal <- ProstateNetworkInputMatrix[which(abs(ProstateNetworkInputMatrix$cor) > 0.95 | ProstateNetworkInputMatrix$p < 0.0000001), c(1,2,3,4)]
write.csv(BreastFinal, "BreastFinalTest.csv")
write.csv(ProstateFinal, "ProstateFinalTest.csv")


##Intersecting transcripts in breast and prostate cancer types as marked in the dataset

BDEG$match <- match(BDEG$location, PDEG$location, nomatch=0)


# Reordering of respective breast cancer and prostate cancer datsets.
# Prostate Normal [1:6], Prostate Tumor [7:12], Breast Normal [13:18], Breast Tumor [19:24]

FinalPDEG <- FinalPDEG [c(4,6,8,10,12,14, 5,7,9,11,13,15, 16,18,20,22,24,26, 17,19,21,23,25,27)]
WorkFinalPDEG <- FinalPDEG[1:12]


FinalBDEG <- FinalBDEG [c(4,6,8,10,12,14, 5,7,9,11,13,15, 16,18,20,22,24,26, 17,19,21,23,25,27)]
WorkFinalBDEG <- FinalBDEG[13:24]

##Prostate data regrerssion analysis(linear model)

par(mfrow=c(1,6))
plot(log2(WorkFinalPDEG$GSM662756),log2(WorkFinalPDEG$GSM662757),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662756) ~ log2(WorkFinalPDEG$GSM662757)), col= 1)
plot(log2(WorkFinalPDEG$GSM662758),log2(WorkFinalPDEG$GSM662759),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662758) ~ log2(WorkFinalPDEG$GSM662759)), col= 1)
plot(log2(WorkFinalPDEG$GSM662760),log2(WorkFinalPDEG$GSM662761),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662760) ~ log2(WorkFinalPDEG$GSM662761)), col= 1)
plot(log2(WorkFinalPDEG$GSM662762),log2(WorkFinalPDEG$GSM662763),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662762) ~ log2(WorkFinalPDEG$GSM662763)), col= 1)
plot(log2(WorkFinalPDEG$GSM662764),log2(WorkFinalPDEG$GSM662765),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662764) ~ log2(WorkFinalPDEG$GSM662765)), col= 1)
plot(log2(WorkFinalPDEG$GSM662766),log2(WorkFinalPDEG$GSM662767),  pch = 16, cex = 1.3, col = c("blue","red"))
abline(lm(log2(WorkFinalPDEG$GSM662766) ~ log2(WorkFinalPDEG$GSM662767)), col= 1)


##Gene Set Enrichment Analysis

library(genefilter)
library(GSEABase)
Breast_GSEA <- GeneSetCollection(WorkFinalBDEG, setType = KEGGCollection())
Prostate_GSEA <- GeneSetCollection(WorkFinalPDEG, setType = KEGGCollection())


# Support Vector Machine Implementation
## Prostate

install.packages("e1071")
library(e1071)
temp1 <- WorkFinalPDEG
temp1$ID_REF <- NULL
temp1 <- log2(temp1)
temp1 <- t(temp1)
ClassLabels1 <- c(rep(1,6),rep(-1,6))
DataFrame1 <- data.frame(Gene=temp1,ClassLabels=as.factor(ClassLabels1))
SVMModel1 <- svm(ClassLabels1~., data=DataFrame1, kernel="linear", cost=10, scale = FALSE)
GeneWeights1<-t(SVMModel1$coefs)%*%SVMModel1$SV
sort.list(GeneWeights1) ## Genes 212 and 129 have highest and second highest weights, respectively. 
plot(SVMModel1,DataFrame1, Gene.212 ~ Gene.129)

##Breast

temp2 <- WorkFinalBDEG
temp2$ID_REF <- NULL
temp2 <- log2(temp2)
temp2 <- t(temp2)
ClassLabels2<- c(rep(1,6),rep(-1,6))
DataFrame2 <- data.frame(Gene=temp2,ClassLabels=as.factor(ClassLabels2))
SVMModel2 <- svm(ClassLabels2~., data=DataFrame2, kernel="linear", cost=10, scale = FALSE)
GeneWeights2<-t(SVMModel2$coefs)%*%SVMModel2$SV
sort.list(GeneWeights2) ## Genes 346 and 133 have highest and second highest weights, respectively. 
plot(SVMModel2,DataFrame2, Gene.346 ~ Gene.133)


install.packages("kernlab")
library(kernlab)
x <- as.matrix(temp1)
y <- matrix(c(rep(1,6),rep(-1,6)))
svp <- ksvm(x,y,type="C-svc", prob.model= TRUE)
predict (svp, x, type= "probabilities")




