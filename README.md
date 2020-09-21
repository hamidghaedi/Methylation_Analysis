# methylation_analysis
Comprehensive tutorial for differential methylation analysis, differential variability analysis and integrative analysis.
The code and approaches that I share here are those I am using to analyze TCGA methylation data.

The three main analysis on methylation data, are covered:

- differential methylation analysis
- differential variability analysis
- integrative analysis

## Workspace preparation
```R
#______________loading packages__________________#
library(TCGAbiolinks)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)

#_________ increasing the memory limit____________#
memory.limit(size = 28000)
```
## Data download and preparation
```R
# DNA methylation aligned to hg19
query_met <- GDCquery(project= "TCGA-BLCA", 
                           data.category = "DNA methylation", 
                           platform = "Illumina Human Methylation 450", 
                           legacy = TRUE)
GDCdownload(query_met)
#putting files togathers
data.met <- GDCprepare(query_met)
#saving the met object
saveRDS(object = data.met,
        file = "data.met.RDS",
        compress = FALSE)
# loading saved session: Once you saved your data, you can load it into your environment: 
data.met = readRDS(file = "data.met.RDS")
# met matrix
met <- as.data.frame(SummarizedExperiment::assay(data.met))
# clinical data
clinical <- data.frame(data.met@colData)

#___________inspectiong methylation data_______________#

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## remove probes with NA
probe.na <- rowSums(is.na(met))

table(probe.na == 0)
#FALSE   TRUE 
#103553 382024 
# chose those has not NA values in rows
probe <- probe.na[probe.na == 0]
met <- met[row.names(met) %in% names(probe), ]

## remove probes that match to chromosome  X and Y 
keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
met <- met[keep, ]
rm(keep) # remove no further needed probes.

## remove SNPs overlapped probe
table(is.na(ann450k$Probe_rs))
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]

# filtre met
met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]

#remove no-further needed dataset
rm(no.snp.probe, probe, probe.na, snp.probe, snp5.probe)

## Removing probes that have been demonstrated to map to multiple places in the genome.
# list adapted from https://www.tandfonline.com/doi/full/10.4161/epi.23470

crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]

# filtre met
met <- met[ -which(row.names(met) %in% crs.reac), ]
bval <- met

## converting beta values to m_values
## m = log2(beta/1-beta)
mval <- t(apply(met, 1, function(x) log2(x/(1-x))))
#______________saving/loading_____________________#
# save data sets
#saveRDS(mval, file = "mval.RDS", compress = FALSE)
#saveRDS (bval, file = "bval.RDS", compress = FALSE)
#mval <- readRDS("mval.RDS")
#bval <- readRDS("bval.RDS")
```
## Differential methylation analysis
### Different methylated CpGs (DMC) or Probe-wise differential methylation analysis
Here we are intrested to find probes that are diffrentially methylated between low-grade and high-grade bladder cancer. In the ```clinical``` dataset , there is column ```paper_Histologic.grade``` that provides data on tumor grade. Also as one can see, most of low-grade bladder cancer samples belong to the *luminal papillary* subtype. 

```R
table(clinical$paper_Histologic.grade, clinical$paper_mRNA.cluster)

#             Basal_squamous Luminal Luminal_infiltrated Luminal_papillary  ND Neuronal
#  High Grade            143      28                  79               120   4       20
#  Low Grade               1       0                   0                20   0        0
#  ND                      0       0                   1                 2   0        0


# filtering and grouping
clinical <- clinical[, c("barcode", "paper_Histologic.grade", "paper_mRNA.cluster")]
clinical <- na.omit(clinical)
clinical <- clinical[-which(clinical$paper_mRNA.cluster == "ND"), ]
clinical <- clinical[-which(clinical$paper_Histologic.grade == "ND"), ]
clinical <- clinical[which(clinical$paper_mRNA.cluster == "Luminal_papillary"), ]
barcode <- clinical$barcode

# removing samples from meth matrixes
bval <- bval[, colnames(bval) %in% barcode]
mval <- mval[, colnames(mval) %in% barcode]

# Making sure about samples in clinical and matrixes and their order
table(colnames(mval) %in% row.names(clinical))
table(colnames(bval) %in% row.names(clinical))
#
all(row.names(clinical) == colnames(bval))
all(row.names(clinical) == colnames(mval))

#Making grouping variable
clinical$paper_Histologic.grade <- as.factor(clinical$paper_Histologic.grade)
#levels(clinical$paper_Histologic.grade)
clinical$paper_Histologic.grade <- relevel(clinical$paper_Histologic.grade, ref = "Low Grade")

#_____________ DMC analysis________________#
design <- model.matrix(~ paper_Histologic.grade, data = clinical)
# fit the linear model 
fit <- lmFit(mval, design)
fit2 <- eBayes(fit)

# extracting significantly methylated probes
deff.meth = topTable(fit2, coef=ncol(design), sort.by="p",number = nrow(mval), adjust.method = "BY")
```

|           |  logFC  |AveExpr  |    t    |   P.Value   | adj.P.Val   |     B   |
|---------- |:-------:|:-------:|--------:|------------:|------------:|--------:|
|cg16177647 |1.629042 |2.901692 |9.204956 |4.442382e-16 |1.511929e-09 |25.51849 |
|cg15465367 |3.319291 |4.106852 |9.129078 |6.895310e-16 |1.511929e-09 |25.10705 |
|cg04117396 |3.699929 |2.056517 |9.019180 |1.301506e-15 |1.902534e-09 |24.51242 |
|cg23797032 |3.051319 |1.563557 |8.958209 |1.849956e-15 |2.028192e-09 |24.18321 |
|cg17278466 |4.550559 |2.647423 |8.881866 |2.870870e-15 |2.517973e-09 |23.77173 |
|cg10113582 |2.614498 |5.553400 |8.576775 |1.645866e-14 |1.202959e-08 |22.13599 |

```R
# Visualization
# plot the top 10 most significantly differentially methylated CpGs 
par(mfrow=c(2,5))
sapply(rownames(deff.meth)[1:10], function(cpg){
  plotCpg(bval, cpg=cpg, pheno=clinical$paper_Histologic.grade, ylab = "Beta values")
})
```
![alt text](https://github.com/hamidghaedi/methylation_analysis/blob/master/dmc.png) 

```R
# making a volcano plot
#making dataset
dat <- data.frame(foldchange = fit[["coefficients"]][,2], logPvalue =  -log10(fit2[["p.value"]][,2]))
dat$threshold <- as.factor(abs(dat$foldchange) < 0.4)

#Visualization
ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
  geom_point(alpha=.6, size=1.2) +
  theme(legend.position="none") +
  xlab("Fold Change") +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none")
  ```
  ![alt text](https://github.com/hamidghaedi/methylation_analysis/blob/master/volcano.plot.PNG)
  
  
### Differentially methylated regions (DMRs)

