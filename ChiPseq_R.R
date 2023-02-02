# CHIPseq analysis
#The BiocManager package, as the modern successor package to BiocInstaller, 
#allows users to install and manage packages from the Bioconductor project. 
#Bioconductor focuses on the statistical analysis and comprehension of high-throughput genomic data.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DiffBind")
BiocManager::install("ChIPseeker")
BiocManager::install("GenomicRanges")
BiocManager::install("ReactomePA")
BiocManager::install("clusterProfiler")
BiocManager::install("TxDb.Hsapiens.UCSC.hg18.knownGene")
BiocManager::install("org.Hs.eg.db")
#BiocManager::install(c("GenomicRanges", "SummarizedExperiment")) - way to install multiple packages at once


##part 1: load the packages in the workspace
library(DiffBind)
library(ChIPseeker)
library(GenomicRanges)
library(ReactomePA)
library(clusterProfiler)

##part 2: import the data into the RStudio workspace
tamoxifen<-dba(sampleSheet="tamoxifen.csv", dir=system.file("extra",package="DiffBind"))
#tamoxifen<-dba(sampleSheet="tamoxifen.csv", dir="PATH/TO/SPREADSHEET") #### to import your own data
##part 3: differential analysis of peaks
######## take a first look at the data
plot(tamoxifen)

######## obtain the counts
tamoxifen.counts <- dba.count(tamoxifen)

######## look at the correlation among the samples after reading counts
plot(tamoxifen.counts)

######## perform differential analysis to look at binding loci enriched in resistant vs. responsive lines
contrast<-dba.contrast(tamoxifen.counts, reorderMeta=list(Condition="Responsive"))#we put the control, here responsive cells
#contrast <- dba.analyze(contrast)
#Blacklist is the list of genes that we don't use for the analysis 
contrast <- dba.analyze(contrast, bGreylist = F, bBlacklist = F) ## if error with greylist/blacklist appears

tamoxifen.DB<-dba.report(contrast)

######## check the number of loci with logFoldChange>0
sum(tamoxifen.DB$Fold>0)
######## check the number of loci with logFoldChange<0
sum(tamoxifen.DB$Fold<0)




##Part 4: visualize the data
######## Principal Component Analysis (PCA)
dba.plotPCA(contrast,DBA_TISSUE,label=DBA_CONDITION)
dba.plotPCA(contrast, contrast=1, label=DBA_TISSUE) ##to see the separation of the diff. binding sites

######## Volcano plot
dba.plotVolcano(contrast)



##Part 5: Annotation and Gene Enrichment of the data
######## conversion of the data to GRanges using GenomicRanges
data = as.data.frame(tamoxifen.DB)
data.up = data[which(data$Fold>0),]
data.down = data[which(data$Fold<0),]
data.up <- GRanges(seqnames=data.up[,1],
             ranges=IRanges(start=data.up[,2],
                            end= data.up[,3],
                            names=paste("Site", 1:nrow(data.up), sep="")),
             Fold = data.up$Fold,
             FDR = data.up$FDR)
data.down <- GRanges(seqnames=data.down[,1],
                   ranges=IRanges(start=data.down[,2],
                                  end= data.down[,3],
                                  names=paste("Site", 1:nrow(data.down), sep="")),
                   Fold = data.down$Fold,
                   FDR = data.down$FDR)

#X = GRanges(X)
######## Annotation of the data using ChIPSeeker
#BiocManager::install("TxDb.Hsapiens.UCSC.hg18.knownGene")
#BiocManager::install("org.Hs.eg.db")
txdb = TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene
anno_data.up <- as.data.frame(annotatePeak(data.up, tssRegion=c(-3000, 3000),
                                                     TxDb=txdb, annoDb="org.Hs.eg.db"))
anno_data.down <- as.data.frame(annotatePeak(data.down, tssRegion=c(-3000, 3000),
                                           TxDb=txdb, annoDb="org.Hs.eg.db"))
#OR
peakAnnoList <- lapply(list(data.up, data.down), annotatePeak, TxDb=txdb,
                       tssRegion=c(-200, 200), verbose=FALSE)

######## Visualization of binding features
plotAnnoBar(peakAnnoList)

######## Gene Ontology (GO) analysis
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
data = data.frame("Entrez" = c(genes[[1]], genes[[2]]),"group" = c(rep(1,length(genes[[1]])), rep(2,length(genes[[2]]))))
compGO <- compareCluster(geneCluster = Entrez~group,
                         fun = "enrichGO",
                         data = data,
                         OrgDb = 'org.Hs.eg.db', pvalueCutoff=0.1, ont = "BP")
dotplot(compGO, title = "GO Analysis")
