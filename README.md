# Methylation Analysis
Comprehensive tutorial for differential methylation analysis, differential variability analysis and integrative analysis.
The code and approaches that I share here are those I am using to analyze TCGA methylation data.

The three main analysis on methylation data, are covered:

- differential methylation analysis
- differential variability analysis
- integrative analysis *to be aaded

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
library(edgeR)

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
  
  
### Differentially methylated regions (DMRs) analysis
```R
# setting some annotation
myAnnotation <- cpg.annotate(object = mval, datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = "paper_Histologic.gradeHigh Grade", 
                             arraytype = "450K",
                             fdr = 0.001)
str(myAnnotation)

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
```
 ![alt text](https://github.com/hamidghaedi/Methylation_Analysis/blob/master/top.20.dmr.PNG) 

```R
# visualization
dmr.table <- data.frame(results.ranges)

# setting up variable for groupinh and color

pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(clinical$paper_Histologic.grade))]
names(groups) <- levels(factor(clinical$paper_Histologic.grade))

#setting up the genomic region 
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 2
# coordinates are stored under results.ranges[dmrIndex]

chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))

# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))


# defining CpG islands track
# download cpgislands for chromosome number 6 from ucsc
chr6.cpg <- read.csv("chr6-cpg.csv")

islandData <- GRanges(seqnames=Rle(chr6.cpg[,1]), 
                      ranges=IRanges(start=chr6.cpg[,2],
                                     end=chr6.cpg[,3]),
                      strand=Rle(strand(rep("*",nrow(chr6.cpg)))))

# DNAseI hypersensitive sites track
#downloaded from ucsc
chr6.dnase <- read.csv("chr6-dnase.csv")

dnaseData <- GRanges(seqnames=chr6.dnase[,1],
                     ranges=IRanges(start=chr6.dnase[,2], end=chr6.dnase[,3]),
                     strand=Rle(rep("*",nrow(chr6.dnase))),
                     data=chr6.dnase[,5])

#Setting up the ideogram, genome and RefSeq tracks 

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=paste0(chrom))
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)

#Ensure that the methylation data is ordered by chromosome and base position.

ann450kOrd <- ann450k[order(ann450k$chr,ann450k$pos),]
bvalOrd <- bval[match(ann450kOrd$Name,rownames(bval)),]

#Create the data tracks:
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bvalOrd)

# methylation data track
methTrack <- DataTrack(range=cpgData, 
                       groups=clinical$paper_Histologic.grade, # change this if your groups are diffrent
                       genome = gen,
                       chromosome=chrom,
                       ylim=c(-0.05,1.05),
                       col=pal,
                       type=c("a","p"), 
                       name="DNA Meth.\n(beta value)",
                       background.panel="white", 
                       legend=TRUE, 
                       cex.title=0.8,
                       cex.axis=0.8, 
                       cex.legend=0.8)

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")


# Set up the track list and indicate the relative sizes of the different tracks. 
# Finally, draw the plot using the plotTracks function
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks

tiff( filename = "dmr.tiff", width = 15, height = 10, units = "in", res = 400)
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
dev.off()

```
  ![alt text](https://github.com/hamidghaedi/methylation_analysis/blob/master/dmr.png)

## Differential variability analysis

Rather than testing for DMCs and DMRs, we may be interested in testing for differences between group variances. This could be quite helpful for feature selection in ML baesd projects. In these situation you may perefre to selct variables that shows great diffrences between groups.
```R
#__________________________Differential variability_________________#
fitvar <- varFit(mval, design = design, coef = c(1,2))

# Summary of differential variability
summary(decideTests(fitvar))

topDV <- topVar(fitvar)
# Top 10 differentially variable CpGs between old vs. newborns
topDV

# visualization
# get beta values for ageing data
par(mfrow=c(5,2))
sapply(rownames(topDV)[1:10], function(cpg){
  plotCpg(bval, cpg=cpg, pheno= clinical$paper_Histologic.grade, 
          ylab = "Beta values")
})


## if you got this error: Error in plot.new() : figure margins too large 
#Do the following and try again:
#graphics.off()
#par("mar")
#par(mar=c(1,1,1,1))
```
  ![alt text](https://github.com/hamidghaedi/methylation_analysis/blob/master/VMCs.PNG)

## Integrated analysis
Aberrant methylation of CpGs in promoter  potentially could alter downstream gene expression level by interfering transcription factor function, and this process is knows as *cis-regulation*. In the other hand,  CpGs may also affect expression of genes that are located far away from them in genome. Usually this could be observed when if CpGs position  locations are happened to be in regulatory elements like enhancers. This type of regulation between CpG and gene is known as *trans-regulation*.
Steps toward cis-regulation analysis:
- Gene expression analysis (making normalized count matrix, finding diffrentially exxpressed genes)
- finding probes resided within promoters
- finding diffrentially methylated probes (|beta value| > 0.2)
- Performing correlation analysis to find those with |r| > 0.3 and p-value < 0.05

```R
###gene expression data download and analysis
## expression data
query.exp <- GDCquery(project = "TCGA-BLCA",
                      platform = "Illumina HiSeq",
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification", 
                      file.type = "results",
                      legacy = TRUE)
GDCdownload(query.exp,  method = "api")
dat<- GDCprepare(query = query.exp, save = TRUE, save.filename = "blcaExp.rda")

rna <-assay(dat)
clinical.exp = data.frame(colData(dat))

# find what we have for grade
table(clinical.exp$paper_Histologic.grade)
#High Grade  Low Grade         ND 
#384         21                3 
table(clinical.exp$paper_Histologic.grade, clinical.exp$paper_mRNA.cluster)

# Get rid of ND and NA samples, normal samples

clinical.exp <- clinical.exp[(clinical.exp$paper_Histologic.grade == "High Grade" | 
                        clinical.exp$paper_Histologic.grade == "Low Grade"), ]
clinical.exp$paper_Histologic.grade[clinical.exp$paper_Histologic.grade == "High Grade"] <- "High_Grade"
clinical.exp$paper_Histologic.grade[clinical.exp$paper_Histologic.grade == "Low Grade"] <- "Low_Grade"
# since most of low-graded are in Luminal_papilary category, we remain focus on this type
clinical.exp <- clinical.exp[clinical.exp$paper_mRNA.cluster == "Luminal_papillary", ]
clinical.exp <- clinical.exp[!is.na(clinical.exp$paper_Histologic.grade), ]


# keep samples matched between clinical.exp file and expression matrix
rna <- rna[, row.names(clinical.exp)]
all(rownames(clinical.exp) %in% colnames(rna))
#TRUE

## A pipeline for normalization and gene expression analysis (edgeR and limma)

edgeR_limma.pipe = function(
  exp_mat,
  groups,
  ref.group=NULL){  
  group = factor(clinical.exp[, groups])
  if(!is.null(ref.group)){group = relevel(group, ref=ref.group)}
  # making DGEList object
  d = DGEList(counts= exp_mat,
              samples=clinical.exp,
              genes=data.frame(rownames(exp_mat)))
  
  # filtering
  keep = filterByExpr(d,design)
  d = d[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Calculate normalization factors to scale the raw library sizes (TMM and voom)
  design = model.matrix(~ group)
  d = calcNormFactors(d, method="TMM")
  v = voom(d, design, plot=TRUE)
  
  # Model fitting and DE calculation 
  fit = lmFit(v, design)
  fit = eBayes(fit)
  # DE genes
  DE = topTable(fit, coef=ncol(design), sort.by="p",number = nrow(rna), adjust.method = "BY")
  return(
    list( 
      DE=DE, # DEgenes
      voomObj=v, # NOrmalized counts
      fit=fit # DE stats
    )
  )
}

# Runing the pipe
de.list <- edgeR_limma.pipe(rna,"paper_Histologic.grade", "Low_Grade" )
de.genes <- de.list$DE
#ordering diffrentially expressed genes
de.genes<-de.genes[with(de.genes, order(abs(logFC), adj.P.Val, decreasing = TRUE)), ]
# voomObj is normalized expression values on the log2 scale
norm.count <- data.frame(de.list$voomObj)
norm.count <- norm.count[,-1]
norm.count <- t(apply(norm.count,1, function(x){2^x}))
colnames(norm.count) <- chartr(".", "-", colnames(norm.count))

#______________preparing methylation data for cis-regulatory analysis____________#

# finding probes in promoter of genes
table(data.frame(ann450k)$Regulatory_Feature_Group) ## to find regulatory features of probes

# selecting a subset of probes associated with promoted
promoter.probe <- rownames(data.frame(ann450k))[data.frame(ann450k)$Regulatory_Feature_Group 
                                                 %in% c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")]

# find genes probes with significantly different methylation status in 
# low- and high-grade bladder cancer

low.g_id <- clinical$barcode[clinical$paper_Histologic.grade == "Low Grade"]
high.g_id <- clinical$barcode[clinical$paper_Histologic.grade == "High Grade"]

dbet <- data.frame (low.grade = rowMeans(bval[, low.g_id]),
                     high.grade = rowMeans(bval[, high.g_id]))
dbet$delta <- abs(dbet$low.grade - dbet$high.grade)

db.probe <- rownames(dbet)[dbet$delta > 0.2] # those with deltabeta > 0.2
db.probe <- db.probe %in% promoter.probe # those resided in promoter
rm(dbet)

# those genes flanked to promote probe
db.genes <- data.frame(ann450k)[rownames(data.frame(ann450k)) %in% db.probe, ]
db.genes <- db.genes[, c("Name","UCSC_RefGene_Name")]
db.genes <- tidyr::separate_rows(db.genes, Name, UCSC_RefGene_Name) # extending collapsed cells
db.genes$comb <- paste(db.genes$Name,db.genes$UCSC_RefGene_Name) # remove duplicates
db.genes <- db.genes[!duplicated(db.genes$comb), ]
db.genes <- db.genes[, -3]

# doing correlation analysis
# polishing matrices to have only high grade samples
cis.bval.mat <- bval[, high.g_id]
cis.exp.mat <- norm.count[, high.grade.exp.id]
#making patient name similar
colnames(cis.bval.mat) <- substr(colnames(cis.bval.mat),1,19)
colnames(cis.exp.mat) <- substr(colnames(cis.exp.mat),1,19)
cis.exp.mat <- cis.exp.mat[, colnames(cis.bval.mat)]

#editing expression matrix rowname
df <- data.frame(name = row.names(cis.exp.mat)) # keeping rownames as a temporary data frame
df <- data.frame(do.call('rbind', strsplit(as.character(df$name),'|',fixed=TRUE))) # this do magic like "text to column" in Excel!
rowName <- df$X1
# find duplicates in rowName, if any
table(duplicated(rowName))
#FALSE  TRUE 
#20530     1 
# in order to resolve  duplucation issue
rowName[duplicated(rowName) == TRUE]
#[1] "SLC35E2"
#
rowName[grep("SLC35E2", rowName)[2]] <- "SLC35E2_2"
#setting rna row names 
row.names(cis.exp.mat) <- rowName
rm(df, rowName) # removing datasets that we do not need anymore

#__________________correlation analysis__________________#
cis.reg = data.frame( gene=character(0), cpg=character(0), pval=numeric(0), cor=numeric(0))

for (i in 1:nrow(db.genes)){
  cpg = db.genes[i,][1]
  gene = db.genes[i,][2]
  if (gene %in% rownames(cis.exp.mat)){
    df1 <- data.frame(exp= cis.exp.mat[as.character(gene), ])
    df2 <- t(cis.bval.mat[as.character(cpg), ])
    df <- merge(df1,df2, by = 0)
    res <- cor.test(df[,2], df[,3], method = "pearson")
    pval = round(res$p.value, 4)
    cor = round(res$estimate, 4)
    cis.reg[i,] <- c(gene, cpg, pval, cor)
  }
}
cis.reg$adj.P.Val = round(p.adjust(cis.reg$pval, "fdr"),4)
cis.reg <- cis.reg[with(cis.reg, order(cor, adj.P.Val)), ]

# top pair visualization
# inspecting the results, C2orf74 gene has significant correlation with probes:

gen.vis <- merge(data.frame(exp= cis.exp.mat["C2orf74", ]), 
            t(cis.bval.mat[c("cg24757310", "cg01648237", "cg05037927", "cg16328106", "cg23405039", "cg18158151"), ]),
            by = 0)

par(mfrow=c(3,2))
sapply(names(gen.vis)[3:8], function(cpg){
  plot(x= gen.vis[ ,cpg], y = gen.vis[,2], xlab = "beta value", 
       ylab = "normalized expression" ,
       pch = 19,
       main = paste("C2orf74",cpg, sep = "-"),
       frame = FALSE)
  abline(lm(gen.vis[,2] ~ gen.vis[ ,cpg], data = gen.vis), col = "blue")
})
```
![alt text](https://github.com/hamidghaedi/Methylation_Analysis/blob/master/cis_reg.PNG)

