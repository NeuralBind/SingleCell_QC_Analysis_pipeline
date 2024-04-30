#save.image("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/scRNAseq_workflow_RasVCTL_NASH_CANIL.RData")
load("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/scRNAseq_workflow_RasVCTL_NASH_CANIL.RData")

library(dplyr)
library(plyr)
library(Seurat)
library(patchwork)
library(scater)
library(biomaRt)
library(scater)
library(robustbase)
library(cowplot)
library(ggplot2)
library(reshape2)
library(gdata)
library(purrr)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ArchR)
set.seed(23)
############################ 1. Open the datasets  ############################ 
############################ 
OLDCtl <- Read10X_h5(filename = "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_cellranger/feature_bc_matrix/765960_CTR_EV_12_5ug_sample_feature_bc_matrix.h5", use.names = T, unique.features = T)

Ctl_1 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7332_S1_GEX")
Ctl_2 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7332_S2_GEX")

Intermediate_1 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7329_S1_GEX")
Intermediate_2 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7329_S2_GEX")
Intermediate_3 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7329_S3_GEX")

Endpoint_1 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7257_S1_GEX")
Endpoint_2 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7265_S1_GEX")
Endpoint_3 <- Read10X(data.dir =   "/Volumes/Group Akkari/Individual folders/Masami/HCC/scRNAseq_RasVCANIL/7273_S1_GEX")


############################ 

############################# 2. Initialize the Seurat object with the raw (non-normalized data).############################ 
############################ 
S1_OLDCtl_so <- CreateSeuratObject(counts = OLDCtl, project = "Ctl", min.cells = 3, min.features = 200)
S1_OLDCtl_so#8393
S1_OLDCtl_so$replicate <-"R1"
S1_OLDCtl_so$stage <-"Intermediate"
remove(OLDCtl)


S1_Ctl_so <- CreateSeuratObject(counts = Ctl_1, project = "Ctl", min.cells = 3, min.features = 200)
S1_Ctl_so#14167
S1_Ctl_so$replicate <-"R1"
S1_Ctl_so$stage <-"Intermediate"
remove(Ctl_1)

S2_Ctl_so <- CreateSeuratObject(counts = Ctl_2, project = "Ctl", min.cells = 3, min.features = 200)
S2_Ctl_so#11180
S2_Ctl_so$replicate <-"R2"
S2_Ctl_so$stage <-"Intermediate"
remove(Ctl_2)

S1_Intermediate_so <- CreateSeuratObject(counts = Intermediate_1, project = "RasV", min.cells = 3, min.features = 200)
S1_Intermediate_so#8401
S1_Intermediate_so$replicate <-"R1"
S1_Intermediate_so$stage <-"Intermediate"
remove(Intermediate_1)

S2_Intermediate_so <- CreateSeuratObject(counts = Intermediate_2, project = "RasV", min.cells = 3, min.features = 200)
S2_Intermediate_so#12159
S2_Intermediate_so$replicate <-"R2"
S2_Intermediate_so$stage <-"Intermediate"
remove(Intermediate_2)

S3_Intermediate_so <- CreateSeuratObject(counts = Intermediate_3, project = "RasV", min.cells = 3, min.features = 200)
S3_Intermediate_so#7451
S3_Intermediate_so$replicate <-"R3"
S3_Intermediate_so$stage <-"Intermediate"
remove(Intermediate_3)


S1_Endpoint_so <- CreateSeuratObject(counts = Endpoint_1, project = "RasV", min.cells = 3, min.features = 200)
S1_Endpoint_so#9954
S1_Endpoint_so$replicate <-"R1"
S1_Endpoint_so$stage <-"Endpoint"
remove(Endpoint_1)

S2_Endpoint_so <- CreateSeuratObject(counts = Endpoint_2, project = "RasV", min.cells = 3, min.features = 200)
S2_Endpoint_so#10124
S2_Endpoint_so$replicate <-"R2"
S2_Endpoint_so$stage <-"Endpoint"
remove(Endpoint_2)

S3_Endpoint_so <- CreateSeuratObject(counts = Endpoint_3, project = "RasV", min.cells = 3, min.features = 200)
S3_Endpoint_so#10509
S3_Endpoint_so$replicate <-"R3"
S3_Endpoint_so$stage <-"Endpoint"
remove(Endpoint_3)



############################# 

############################## 3. QC and selecting cells for further analysis #############################
############################## 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#Calculate mitochondrial percent

S1_OLDCtl_so[["percent.mt"]] <- PercentageFeatureSet(S1_OLDCtl_so, pattern = "^mt-")
S1_Ctl_so[["percent.mt"]] <- PercentageFeatureSet(S1_Ctl_so, pattern = "^mt-")
S2_Ctl_so[["percent.mt"]] <- PercentageFeatureSet(S2_Ctl_so, pattern = "^mt-")

S1_Intermediate_so[["percent.mt"]] <- PercentageFeatureSet(S1_Intermediate_so, pattern = "^mt-")
S2_Intermediate_so[["percent.mt"]] <- PercentageFeatureSet(S2_Intermediate_so, pattern = "^mt-")
S3_Intermediate_so[["percent.mt"]] <- PercentageFeatureSet(S3_Intermediate_so, pattern = "^mt-")

S1_Endpoint_so[["percent.mt"]] <- PercentageFeatureSet(S1_Endpoint_so, pattern = "^mt-")
S2_Endpoint_so[["percent.mt"]] <- PercentageFeatureSet(S2_Endpoint_so, pattern = "^mt-")
S3_Endpoint_so[["percent.mt"]] <- PercentageFeatureSet(S3_Endpoint_so, pattern = "^mt-")

#View QC metrics
#p1<-VlnPlot(S1_Ctl_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p2<-VlnPlot(MycP53_1_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p3<-VlnPlot(MycP53_2_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p4<-VlnPlot(MycP53_3_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p5<-VlnPlot(MycPten_1_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p6<-VlnPlot(MycPten_2_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p7<-VlnPlot(MycPten_3_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p8<-VlnPlot(Nras12D_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p9<-VlnPlot(Nras12V_1_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p10<-VlnPlot(Nras12V_2_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)
#p11<-VlnPlot(Nras12V_3_so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, raster = F)

#ggpubr::ggarrange(p1,p2,p3,p4,ncol=1,nrow = 4, labels = c("ctl","MycP53_1","MycP53_2","MycP53_3"))
#ggpubr::ggarrange(p5,p6,p7,p8, ncol=1,nrow = 4,  labels = c("MycPten_1","MycPten_2","MycPten_3","Nras12D"))
#ggpubr::ggarrange(p9,p10,p11,ncol=1,nrow = 4,  labels = c("Nras12V_1","Nras12V_2","Nras12V_3"))

#nFeature_RNA is the number of genes detected in each cell
#nCount_RNA is the total number of molecules detected within a cell.


#Combine all objects into one.
RasVCtl_stages.list <- c(S1_Ctl_so, S2_Ctl_so, S1_OLDCtl_so,S1_Intermediate_so, S2_Intermediate_so,S3_Intermediate_so,
                                 S1_Endpoint_so,S2_Endpoint_so,S3_Endpoint_so)
names(RasVCtl_stages.list) <- c("Ctl_1","Ctl_2","OLDCtl_1","Intermediate_1","Intermediate_2","Intermediate_3","Endpoint_1","Endpoint_2","Endpoint_3")
#remove(S1_Ctl_so, S2_Ctl_so,
#         S1_Intermediate_so, S2_Intermediate_so,S3_Intermediate_so,
#         S1_Endpoint_so,S2_Endpoint_so,S3_Endpoint_so)

#One could use hard threshold for the library size, number of genes detected and mitochondrial content based on the distributions seen above. These would need vary across runs and the decision making process is somewhat arbitrary. It may therefore be preferable to rely on outlier detection to identify cells that markedly differ from most cells.
#We saw above that the distribution of the QC metrics is close to Normal. Hence, we can detect outliers using the median and the median absolute deviation (MAD) from the median (not the mean and the standard deviation which both are sensitive to outliers).
#For a given metric, an outlier value is one that lies over some number of MADs away from the median. A cell will be excluded if it is an outlier in the part of the range to avoid, for example low gene counts, or high mitochondrial content. For a normal distribution, a threshold defined with a distance of 3 MADs from the median retains about 99% of values.
#The scater function isOutlier can be used to detect outlier cells based on any metric in the colData table. It returns a boolean vector that identifies outliers. By default it will mark any cell that is 3 MADS in either direction from the median as an outlier.


#####Subset/filter cell that dont pass the QC



#QC independently
for(i in names(RasVCtl_stages.list)){
  #Mt percentage
  high_Mito_percent <- isOutlier(RasVCtl_stages.list[[i]]$percent.mt, log = T, type="higher")
  table(high_Mito_percent)
  attr(high_Mito_percent, "thresholds")[2]
  RasVCtl_stages.list[[i]]$high_Mt_per <- high_Mito_percent
  #nFeature_RNA
  both_nFeature_RNA <- isOutlier(RasVCtl_stages.list[[i]]$nFeature_RNA, log = T, type="both")
  table(both_nFeature_RNA)
  attr(both_nFeature_RNA, "thresholds")
  RasVCtl_stages.list[[i]]$both_nFeature_RNA <- both_nFeature_RNA
  
  RasVCtl_stages.list[[i]]<-subset(RasVCtl_stages.list[[i]], subset = nFeature_RNA > attr(both_nFeature_RNA, "thresholds")[1] & nFeature_RNA < attr(both_nFeature_RNA, "thresholds")[2] & percent.mt < attr(high_Mito_percent, "thresholds")[2])
  
}
#How many cells are left?
unlist(lapply(RasVCtl_stages.list, function(x){dim(x@assays$RNA)[2]}))
#Ctl_1          Ctl_2       OLDCtl_1 Intermediate_1 Intermediate_2 Intermediate_3     Endpoint_1     Endpoint_2     Endpoint_3 
#13122          10311           7834           7598          11103           6629           8308           9573          10007 

###########Preprocess for DoubletFinder:


#perform in a loop and also include the Cell cycle scoring:
options(future.globals.maxSize = 15000 * 1024^2)
for (i in names(RasVCtl_stages.list)) {
  #RasVCtl_stages.list[[i]] <- NormalizeData(RasVCtl_stages.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  #RasVCtl_stages.list[[i]] <- CellCycleScoring(RasVCtl_stages.list[[i]], g2m.features=g2m.genes_mouse$MGI.symbol, s.features=s.genes_mouse$MGI.symbol)
  RasVCtl_stages.list[[i]] <- SCTransform(RasVCtl_stages.list[[i]],vars.to.regress = "percent.mt", assay = "RNA", verbose = T)
  RasVCtl_stages.list[[i]] <- RunPCA(RasVCtl_stages.list[[i]])
  RasVCtl_stages.list[[i]] <- RunUMAP(RasVCtl_stages.list[[i]], dims = 1:50)
  RasVCtl_stages.list[[i]] <- FindNeighbors(RasVCtl_stages.list[[i]], reduction = "pca", dims = 1:50)
  RasVCtl_stages.list[[i]] <- FindClusters(RasVCtl_stages.list[[i]], resolution = 0.05)
}

#########Doublet finder

library(DoubletFinder)

#1) Separate the list of samples
S1_Ctl_so_DF<-RasVCtl_stages.list$Ctl_1
S2_Ctl_so_DF<-RasVCtl_stages.list$Ctl_2
S1_OLDCtl_so_DF<-RasVCtl_stages.list$OLDCtl_1


S1_Intermediate_so<-RasVCtl_stages.list$Intermediate_1
S2_Intermediate_so<-RasVCtl_stages.list$Intermediate_2
S3_Intermediate_so<-RasVCtl_stages.list$Intermediate_3


S1_Endpoint_so<-RasVCtl_stages.list$Endpoint_1
S2_Endpoint_so<-RasVCtl_stages.list$Endpoint_2
S3_Endpoint_so<-RasVCtl_stages.list$Endpoint_3



## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#paramSweep_v3
sweep.res.list_S1_Ctl_DF <- paramSweep_v3(S1_Ctl_so_DF, PCs = 1:50, sct = T)
sweep.stats_S1_Ctl_DF <- summarizeSweep(sweep.res.list_S1_Ctl_DF, GT = FALSE)

sweep.res.list_S2_Ctl_DF <- paramSweep_v3(S2_Ctl_so_DF, PCs = 1:50, sct = T)
sweep.stats_S2_Ctl_DF <- summarizeSweep(sweep.res.list_S2_Ctl_DF, GT = FALSE)

sweep.res.list_S1_OLDCtl_DF <- paramSweep_v3(S1_OLDCtl_so_DF, PCs = 1:50, sct = T)
sweep.stats_S1_OLDCtl_DF <- summarizeSweep(sweep.res.list_S1_OLDCtl_DF, GT = FALSE)

sweep.res.list_S1_Intermediate_DF <- paramSweep_v3(S1_Intermediate_so, PCs = 1:50, sct = T)
sweep.stats_S1_Intermediate_DF <- summarizeSweep(sweep.res.list_S1_Intermediate_DF, GT = FALSE)

sweep.res.list_S2_Intermediate_DF <- paramSweep_v3(S2_Intermediate_so, PCs = 1:50, sct = T)
sweep.stats_S2_Intermediate_DF <- summarizeSweep(sweep.res.list_S2_Intermediate_DF, GT = FALSE)

sweep.res.list_S3_Intermediate_DF <- paramSweep_v3(S3_Intermediate_so, PCs = 1:50, sct = T)
sweep.stats_S3_Intermediate_DF <- summarizeSweep(sweep.res.list_S3_Intermediate_DF, GT = FALSE)

sweep.res.list_S1_Endpoint_DF <- paramSweep_v3(S1_Endpoint_so, PCs = 1:50, sct = T)
sweep.stats_S1_Endpoint_DF <- summarizeSweep(sweep.res.list_S1_Endpoint_DF, GT = FALSE)

sweep.res.list_S2_Endpoint_DF <- paramSweep_v3(S2_Endpoint_so, PCs = 1:50, sct = T)
sweep.stats_S2_Endpoint_DF <- summarizeSweep(sweep.res.list_S2_Endpoint_DF, GT = FALSE)

sweep.res.list_S3_Endpoint_DF <- paramSweep_v3(S3_Endpoint_so, PCs = 1:50, sct = T)
sweep.stats_S3_Endpoint_DF <- summarizeSweep(sweep.res.list_S3_Endpoint_DF, GT = FALSE)


#find.pK
bcmvn_S1_Ctl_DF <- find.pK(sweep.stats_S1_Ctl_DF)
bcmvn_S1_Ctl_DF$pK[which.max(bcmvn_S1_Ctl_DF$BCmetric)] #0.29

bcmvn_S2_Ctl_DF <- find.pK(sweep.stats_S2_Ctl_DF)
bcmvn_S2_Ctl_DF$pK[which.max(bcmvn_S2_Ctl_DF$BCmetric)] #0.06

bcmvn_S1_OLDCtl_DF <- find.pK(sweep.stats_S1_OLDCtl_DF)
bcmvn_S1_OLDCtl_DF$pK[which.max(bcmvn_S1_OLDCtl_DF$BCmetric)] #0.3

bcmvn_S1_Intermediate_DF <- find.pK(sweep.stats_S1_Intermediate_DF)
bcmvn_S1_Intermediate_DF$pK[which.max(bcmvn_S1_Intermediate_DF$BCmetric)] #0.17

bcmvn_S2_Intermediate_DF <- find.pK(sweep.stats_S2_Intermediate_DF)
bcmvn_S2_Intermediate_DF$pK[which.max(bcmvn_S2_Intermediate_DF$BCmetric)] #0.3

bcmvn_S3_Intermediate_DF <- find.pK(sweep.stats_S3_Intermediate_DF)
bcmvn_S3_Intermediate_DF$pK[which.max(bcmvn_S3_Intermediate_DF$BCmetric)] #0.24

bcmvn_S1_Endpoint_DF <- find.pK(sweep.stats_S1_Endpoint_DF)
bcmvn_S1_Endpoint_DF$pK[which.max(bcmvn_S1_Endpoint_DF$BCmetric)] #0.17

bcmvn_S2_Endpoint_DF <- find.pK(sweep.stats_S2_Endpoint_DF)
bcmvn_S2_Endpoint_DF$pK[which.max(bcmvn_S2_Endpoint_DF$BCmetric)] #0.13

bcmvn_S3_Endpoint_DF <- find.pK(sweep.stats_S3_Endpoint_DF)
bcmvn_S3_Endpoint_DF$pK[which.max(bcmvn_S3_Endpoint_DF$BCmetric)] #0.23

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
#To estimate the theoretical number of multiplet rate : https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
#1)Check from the sequencing use the cellRanger html: estimated cells
#2) Apply rule of 3 to get estimated multiplet rate e.g. from table 10000 cells have a 8% multiplet rate, thus, 14984 cells will have 11.9% multiplet rate ((14984*0.08)/10000)

#------------
homotypic.prop_S1_Ctl <- modelHomotypic(S1_Ctl_so_DF@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults
homotypic.prop_S2_Ctl <- modelHomotypic(S2_Ctl_so_DF@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults
homotypic.prop_S1_OLDCtl <- modelHomotypic(S1_OLDCtl_so_DF@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults

homotypic.prop_S1_Intermediate <- modelHomotypic(S1_Intermediate_so@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults
homotypic.prop_S2_Intermediate <- modelHomotypic(S2_Intermediate_so@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults
homotypic.prop_S3_Intermediate <- modelHomotypic(S3_Intermediate_so@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults

homotypic.prop_S1_Endpoint <- modelHomotypic(S1_Endpoint_so@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults
homotypic.prop_S2_Endpoint <- modelHomotypic(S2_Endpoint_so@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults
homotypic.prop_S3_Endpoint <- modelHomotypic(S3_Endpoint_so@meta.data$seurat_clusters) ## ex: annotations <- S1_Ctl_so_DF@meta.data$ClusteringResults


nExp_poi_S1_Ctl <- round(0.119*nrow(S1_Ctl_so_DF@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi_S2_Ctl <- round(0.097*nrow(S2_Ctl_so_DF@meta.data))  
nExp_poi_S1_OLDCtl <- round(0.073*nrow(S1_OLDCtl_so_DF@meta.data))  

nExp_poi_S1_Intermediate <- round(0.067*nrow(S1_Intermediate_so@meta.data))  
nExp_poi_S2_Intermediate <- round(0.097*nrow(S2_Intermediate_so@meta.data))  
nExp_poi_S3_Intermediate <- round(0.060*nrow(S3_Intermediate_so@meta.data))  

nExp_poi_S1_Endpoint <- round(0.082*nrow(S1_Endpoint_so@meta.data))  
nExp_poi_S2_Endpoint <- round(0.082*nrow(S2_Endpoint_so@meta.data))  
nExp_poi_S3_Endpoint <- round(0.084*nrow(S3_Endpoint_so@meta.data))  

nExp_poi.adj_S1_Ctl <- round(nExp_poi_S1_Ctl*(1-homotypic.prop_S1_Ctl))
nExp_poi.adj_S2_Ctl <- round(nExp_poi_S2_Ctl*(1-homotypic.prop_S2_Ctl))
nExp_poi.adj_S1_OLDCtl <- round(nExp_poi_S1_OLDCtl*(1-homotypic.prop_S1_OLDCtl))

nExp_poi.adj_S1_Intermediate <- round(nExp_poi_S1_Intermediate*(1-homotypic.prop_S1_Intermediate))
nExp_poi.adj_S2_Intermediate <- round(nExp_poi_S2_Intermediate*(1-homotypic.prop_S2_Intermediate))
nExp_poi.adj_S3_Intermediate <- round(nExp_poi_S3_Intermediate*(1-homotypic.prop_S3_Intermediate))

nExp_poi.adj_S1_Endpoint <- round(nExp_poi_S1_Endpoint*(1-homotypic.prop_S1_Endpoint))
nExp_poi.adj_S2_Endpoint <- round(nExp_poi_S2_Endpoint*(1-homotypic.prop_S2_Endpoint))
nExp_poi.adj_S3_Endpoint <- round(nExp_poi_S3_Endpoint*(1-homotypic.prop_S3_Endpoint))
#------------

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
S1_Ctl_so_DF <- doubletFinder_v3(S1_Ctl_so_DF, PCs = 1:50, pN = 0.25, pK = 0.29, nExp = nExp_poi_S1_Ctl, reuse.pANN = FALSE, sct = T)
S1_Ctl_so_DF <- doubletFinder_v3(S1_Ctl_so_DF, PCs = 1:50, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj_S1_Ctl, reuse.pANN = "pANN_0.25_0.29_1562", sct = T)
S1_Ctl_so_DF$DF<-S1_Ctl_so_DF@meta.data$DF.classifications_0.25_0.16_1686
DimPlot(S1_Ctl_so_DF, reduction = "umap",group.by = "DF.classifications_0.25_0.16_1096", shuffle = T)
table(S1_Ctl_so_DF$DF.classifications_0.25_0.16_1417)
table(S1_Ctl_so_DF$DF.classifications_0.25_0.16_1096)


S2_Ctl_so_DF <- doubletFinder_v3(S2_Ctl_so_DF, PCs = 1:50, pN = 0.25, pK = 0.06, nExp = nExp_poi_S2_Ctl, reuse.pANN = FALSE, sct = T)
S2_Ctl_so_DF <- doubletFinder_v3(S2_Ctl_so_DF, PCs = 1:50, pN = 0.25, pK = 0.06, nExp = nExp_poi.adj_S2_Ctl, reuse.pANN = "pANN_0.25_0.06_1000", sct = T)

S1_OLDCtl_so_DF <- doubletFinder_v3(S1_OLDCtl_so_DF, PCs = 1:50, pN = 0.25, pK = 0.3, nExp = nExp_poi_S1_OLDCtl, reuse.pANN = FALSE, sct = T)
S1_OLDCtl_so_DF <- doubletFinder_v3(S1_OLDCtl_so_DF, PCs = 1:50, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj_S1_OLDCtl, reuse.pANN = "pANN_0.25_0.3_572", sct = T)


S1_Intermediate_so <- doubletFinder_v3(S1_Intermediate_so, PCs = 1:50, pN = 0.25, pK = 0.17, nExp = nExp_poi_S1_Intermediate, reuse.pANN = FALSE, sct = T)
S1_Intermediate_so <- doubletFinder_v3(S1_Intermediate_so, PCs = 1:50, pN = 0.25, pK = 0.17, nExp = nExp_poi.adj_S1_Intermediate, reuse.pANN = "pANN_0.25_0.17_509", sct = T)


S2_Intermediate_so <- doubletFinder_v3(S2_Intermediate_so, PCs = 1:50, pN = 0.25, pK = 0.3, nExp = nExp_poi_S2_Intermediate, reuse.pANN = FALSE, sct = T)
S2_Intermediate_so <- doubletFinder_v3(S2_Intermediate_so, PCs = 1:50, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj_S2_Intermediate, reuse.pANN = "pANN_0.25_0.3_1077", sct = T)


S3_Intermediate_so <- doubletFinder_v3(S3_Intermediate_so, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi_S3_Intermediate, reuse.pANN = FALSE, sct = T)
S3_Intermediate_so <- doubletFinder_v3(S3_Intermediate_so, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj_S3_Intermediate, reuse.pANN = "pANN_0.25_0.24_398", sct = T)


S1_Endpoint_so <- doubletFinder_v3(S1_Endpoint_so, PCs = 1:50, pN = 0.25, pK = 0.17, nExp = nExp_poi_S1_Endpoint, reuse.pANN = FALSE, sct = T)
S1_Endpoint_so <- doubletFinder_v3(S1_Endpoint_so, PCs = 1:50, pN = 0.25, pK = 0.17, nExp = nExp_poi.adj_S1_Endpoint, reuse.pANN = "pANN_0.25_0.17_681", sct = T)


S2_Endpoint_so <- doubletFinder_v3(S2_Endpoint_so, PCs = 1:50, pN = 0.25, pK = 0.13, nExp = nExp_poi_S2_Endpoint, reuse.pANN = FALSE, sct = T)
S2_Endpoint_so <- doubletFinder_v3(S2_Endpoint_so, PCs = 1:50, pN = 0.25, pK = 0.13, nExp = nExp_poi.adj_S2_Endpoint, reuse.pANN = "pANN_0.25_0.13_785", sct = T)


S3_Endpoint_so <- doubletFinder_v3(S3_Endpoint_so, PCs = 1:50, pN = 0.25, pK = 0.23, nExp = nExp_poi_S3_Endpoint, reuse.pANN = FALSE, sct = T)
S3_Endpoint_so <- doubletFinder_v3(S3_Endpoint_so, PCs = 1:50, pN = 0.25, pK = 0.23, nExp = nExp_poi.adj_S3_Endpoint, reuse.pANN = "pANN_0.25_0.23_841", sct = T)


#Combine all objects into one.
RasVCtl_stages.list_DF <- c(S1_Ctl_so_DF, S2_Ctl_so_DF, S1_OLDCtl_so_DF,S1_Intermediate_so, S2_Intermediate_so,S3_Intermediate_so,
                         S1_Endpoint_so,S2_Endpoint_so,S3_Endpoint_so)
names(RasVCtl_stages.list_DF) <- c("Ctl_1","Ctl_2","OLDCtl_1","Intermediate_1","Intermediate_2","Intermediate_3","Endpoint_1","Endpoint_2","Endpoint_3")

#How many doublets were identified:
for(i in names(RasVCtl_stages.list_DF)){
  print(i)
  print(table(RasVCtl_stages.list_DF[[i]]@meta.data[,ncol(RasVCtl_stages.list_DF[[i]]@meta.data)]))
}

#"Ctl_1" -Doublet Singlet 9.06%
#1189   11933 
#"Ctl_2"- Doublet Singlet 7.5%
#775    9536 
#"OLDCtl_1" - Doublet Singlet 30%
#432    7402
#"Intermediate_1"-nDoublet Singlet 5%
#405    7193 
#"Intermediate_2"-Doublet Singlet 8%
#912   10191
#"Intermediate_3"-Doublet Singlet  4%
#316    6313
#"Endpoint_1"-Doublet Singlet  6%
#445    7863 
#"Endpoint_2"-Doublet Singlet 
#597    8976
#"Endpoint_3"-Doublet Singlet 
#650    9357 

#Save list of object with doublets annotated:

save("RasVCtl_stages.list_DF", file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasVCtl_stages.list_withdoublets.RData")
load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasVCtl_stages.list_withdoublets.RData")

#Remove doublets

for(i in names(RasVCtl_stages.list_DF)){
  #subset only singles
  cellids<-rownames(RasVCtl_stages.list_DF[[i]]@meta.data[RasVCtl_stages.list_DF[[i]]@meta.data[,ncol(RasVCtl_stages.list_DF[[i]]@meta.data)]=="Singlet",])
  RasVCtl_stages.list_DF[[i]]<-subset(RasVCtl_stages.list_DF[[i]], cells= cellids)
  
}
RasVCtl_stages.list_DF$OLDCtl_1$stage<-"OLD_Intermediate"
#How many cells are left?
unlist(lapply(RasVCtl_stages.list_DF, function(x){dim(x@assays$RNA)[2]}))
#Ctl_1          Ctl_2       OLDCtl_1 Intermediate_1 Intermediate_2 Intermediate_3     Endpoint_1     Endpoint_2     Endpoint_3 
#11933           9536           7402           7193          10191           6313           7863           8976           9357 

##############################


###############################4. Normalisation:############################## 
############################## 

#Cell Cycle scoring:
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#Transform to mouse genes

#s.genes_mouse<-getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = s.genes, mart = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("hgnc_symbol"), martL = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
#g2m.genes_mouse<-getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = g2m.genes, mart = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("hgnc_symbol"), martL = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)

#perform in a loop and also include the Cell cycle scoring:
options(future.globals.maxSize = 28000 * 1024^2)
for (i in names(RasVCtl_stages.list_DF)) {
  RasVCtl_stages.list_DF[[i]] <- NormalizeData(RasVCtl_stages.list_DF[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  #RasVCtl_stages.list[[i]] <- CellCycleScoring(RasVCtl_stages.list[[i]], g2m.features=g2m.genes_mouse$MGI.symbol, s.features=s.genes_mouse$MGI.symbol)
  RasVCtl_stages.list_DF[[i]] <- SCTransform(RasVCtl_stages.list_DF[[i]],vars.to.regress = "percent.mt", assay = "RNA", verbose = T)
}

############################## 

################################ 5. Identification of highly variable features (feature selection):############################## 
############################## 
for (i in names(RasVCtl_stages.list_DF)) {
  RasVCtl_stages.list_DF[[i]] <- FindVariableFeatures(RasVCtl_stages.list_DF[[i]],  selection.method = "vst", nfeatures = 2000)
}
############################## 




################################ #6. Integrate datasets:################################ 
############################## 
#listgenopheno<-list(ctl_so, Nras12D_so)
#features <- SelectIntegrationFeatures(object.list = listgenopheno)
#listgenopheno <- PrepSCTIntegration(object.list = listgenopheno, anchor.features = features)
#anchors <- FindIntegrationAnchors(object.list = listgenopheno,normalization.method = "SCT", anchor.features = features)

#IMPORTANT UPDATE: OLD CONTROL REMOVED
RasVCtl_stages.list_DF[["OLDCtl_1"]]<-NULL

features <- SelectIntegrationFeatures(object.list = RasVCtl_stages.list_DF)
save(features, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasV_Ctl_stages_features.RData")
#load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasV_Ctl_stages_features.RData")

RasVCtl_stages.list_DF <- PrepSCTIntegration(object.list = RasVCtl_stages.list_DF, anchor.features = features)
for (i in names(RasVCtl_stages.list_DF)) {
  #RasVCtl_stages.list[[i]] <- ScaleData(RasVCtl_stages.list[[i]],  features = features, verbose = FALSE)
  RasVCtl_stages.list_DF[[i]] <- RunPCA(RasVCtl_stages.list_DF[[i]],  features = features, verbose = T)
}
ElbowPlot(RasVCtl_stages.list_DF[[i]], ndims = 50)
anchors <- FindIntegrationAnchors(object.list = RasVCtl_stages.list_DF,normalization.method = "SCT", anchor.features = features,reduction = "rpca",
                                  dims = 1:50)


# this command creates an 'integrated' data assay
#combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
RasVCtl_stages.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)
remove(anchors)


############################## 




############################## #7. Dimension reduction and clustering ################################ 
############################## 

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay

DefaultAssay(RasVCtl_stages.combined) <- "integrated"
RasVCtl_stages.combined <- RunPCA(RasVCtl_stages.combined, verbose = T)
#DimPlot(RasVCtl_stages.combined, reduction = "pca")


ElbowPlot(RasVCtl_stages.combined, ndims = 50)
RasVCtl_stages.combined <- RunUMAP(RasVCtl_stages.combined, reduction = "pca", dims = 1:50)
DimPlot(RasVCtl_stages.combined, reduction = "umap", shuffle = T)

RasVCtl_stages.combined <- FindNeighbors(RasVCtl_stages.combined, reduction = "pca", dims = 1:50)
RasVCtl_stages.combined <- FindClusters(RasVCtl_stages.combined, resolution = 0.05)
  

#save(features, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasV_Ctl_stages_features.RData")
#load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasV_Ctl_stages_features.RData")

#save(RasVCtl_stages.combined, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasV_Ctl_stages.RData")
#load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/RasV_Ctl_stages_features.RData")


RasVCtl_stages.combined$sample_name<-paste(RasVCtl_stages.combined$orig.ident,RasVCtl_stages.combined$stage, sep = "_" )
#RasVCtl_stages.combined$sample_name<-gsub(pattern = "Ctl_OLD_Intermediate","Ctl",RasVCtl_stages.combined$sample_name )
RasVCtl_stages.combined$sample_name<-gsub(pattern = "Ctl_Intermediate","Ctl",RasVCtl_stages.combined$sample_name)
RasVCtl_stages.combined$sample_name<-gsub(pattern = "RasV_Intermediate","RasV_Pretumor",RasVCtl_stages.combined$sample_name )

cols<-ArchRPalettes[[8]][c(3,10,11)]
names(cols)<-names(table(RasVCtl_stages.combined$sample_name))
DimPlot(RasVCtl_stages.combined, cols = cols[!is.na(names(cols))], reduction = "umap",label = F, repel = TRUE, shuffle = T, group.by = "sample_name")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_Ras12V_samplename_res05.pdf", width = 6, height = 5)
# Visualization
#p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
#p2 <- DimPlot(combined, reduction = "umap",label = TRUE, repel = TRUE)
cols<-rev(ArchRPalettes[[9]])
cols2<-ArchRPalettes[[11]][c(1,3)]
#cols2<-c("grey","#FF8C00")
names(cols2)<-names(table(RasVCtl_stages.combined$orig.ident))
names(cols)<-names(table(Idents(RasVCtl_stages.combined)))

p1 <- DimPlot(RasVCtl_stages.combined, cols = cols2[!is.na(names(cols2))], reduction = "umap", group.by = "orig.ident", shuffle = T)
p2 <- DimPlot(RasVCtl_stages.combined, reduction = "umap", group.by = "stage", shuffle = T)
p3 <- DimPlot(RasVCtl_stages.combined, cols = cols[!is.na(names(cols))], reduction = "umap",label = F, repel = TRUE, shuffle = T)
p3
p1 + p2 + p3
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_Ras12V_res025.pdf", width = 15, height = 5)
#DimPlot(combined, reduction = "umap", split.by = "orig.ident")

p3
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_Ras12V_res05_annotated.pdf", width = 7, height = 5)

#Remove CD45neg cells:Still has cells in the endothelial cluster
#RasVCtl_stages.combinedminusCD45neg<-subset(x = RasVCtl_stages.combined, subset = Ptprc>0)
#DefaultAssay(RasVCtl_stages.combined) <- "SCT"
FeaturePlot(RasVCtl_stages.combined, features = c("Ptprc")) #Cd45+ =Leukocytes cells; #Cd11b+ =Myeloid cells

############################## 

############################### 8. Plot markers ###############################
############################## 
DefaultAssay(RasVCtl_stages.combined) <- "SCT"
FeaturePlot(RasVCtl_stages.combined, features = c("Ly6c1","Ly6c2","Itgam","Adgre1","Ly6g","Csf2ra","Csf2rb")) #Cd45+ =Leukocytes cells; #Cd11b+ =Myeloid cells
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/HDTV_RasD/Featureplot_Ly6Clow_F480Low_Csf2rab.pdf", height = 12, width = 13)



FeaturePlot(RasVCtl_stages.combined, features = c("Plvap")) #Cd45+ =Leukocytes cells; #Cd11b+ =Myeloid cells

RasVCtl_stages.combined[["orig.ident"]] <- factor(x = RasVCtl_stages.combined@meta.data$orig.ident, levels = c("IgG","GMCSF"))
names(cols2)<-names(table(RasVCtl_stages.combined$orig.ident))
VlnPlot(RasVCtl_stages.combined,fill.by = "ident",cols = cols2[!is.na(names(cols2))],stack = T, flip = T,  split.by  = "orig.ident",features =c("Entpd1","Cd80","Cd40","Cd274","H2-Ab1","Itgax","Arg1","Nos2","Il1a","Il1b","Il6","Ccl17","Ccl6","Csf1r","Cx3cr1","Mrc1","Ccr2"),slot = "data", assay = "RNA", pt.size = 0.2)
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/HDTV_RasD/Vln_Ly6Clow_F480Low_Csf2rab_InflImmuno.pdf", height = 15, width = 13)


DefaultAssay(RasVCtl_stages.combined) <- "SCT"
FeaturePlot(RasVCtl_stages.combined, features = c("Csf2")) #Cd45+ =Leukocytes cells; #Cd11b+ =Myeloid cells
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/HDTV_RasD/Featureplot_Csf2.pdf", height = 5, width = 5)


FeaturePlot(clust0_Monocytic,features ="G2M.Score" ) #Cd45+ =Leukocytes cells; #Cd11b+ =Myeloid cells
FeaturePlot(RasVCtl_stages.combined, features = c("Ly6c2","Ly6g")) 
DefaultAssay(RasVCtl_stages.combinedminusEndothelial) <- "SCT"
FeaturePlot(RasVCtl_stages.combinedminusEndothelial, features = c("Ceacam2")) #C
FeaturePlot(RasVCtl_stages.combinedminusEndothelial, features = c("Siglecf")) #C

VlnPlot(RasVCtl_stages.combined, features =c("Arg1","Csf2ra","Ccl6"),slot = "data", assay = "RNA",group.by = "orig.ident", pt.size = 0.2)
VlnPlot(RasVCtl_stages.combined,fill.by = "ident",cols = cols2[!is.na(names(cols2))],stack = T, flip = T,  split.by  = "orig.ident",features =c("Arg1","Csf2rb","Ccl17","Ccl6","Ms4a3"),slot = "data", assay = "RNA", pt.size = 0.2)
VlnPlot(RasVCtl_stages.combined, features =c("Irf8","Klf4","Nr4a1","Zbtb46"),stack = T, flip = T, slot = "data", assay = "RNA",split.by  = "orig.ident", pt.size = 0.2)
ggsave("/Users/m.ando/surfdrive/Documents/HCC/Figures/Hannahgenes_RasDCtl_Vlnplot.pdf", height = 5, width = 8)

FeaturePlot(RasVCtl_stages.combined, features =c("Irf8","Klf4","Nr4a1","Zbtb46"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/Figures/Hannahgenes_RasDCtl_Featureplot.pdf", height = 6, width = 7)


colorsfeatureplots<-ArchRPalettes$whitePurple
FeaturePlot(RasVCtl_stages.combined, features = c("Il12a","Il12b","Csf2","Csf2ra","Csf2rb","Ccl6"), cols=colorsfeatureplots) #Cd11b+ =Myeloid cells

#GMSCF signalling 
FeaturePlot(RasVCtl_stages.combined, features = c("Csf2ra","Csf2rb","Ccl6","Ccl24","Ccl17")) #Cd11b+ =Myeloid cells
FeaturePlot(RasVCtl_stages.combined, features = c("Il2","Il12a","Nt5e","Trem2")) #Cd11b+ =Myeloid cells

#Change order of pops
Idents(RasVCtl_stages.combined) <- factor(x = Idents(RasVCtl_stages.combined), levels = levels(RasVCtl_stages.combined)[c(1,3,2,4:8)])
#Pro Inflammatory
VlnPlot(RasVCtl_stages.combined,idents =levels(RasVCtl_stages.combined)[c(1,3,2,4:5,7:8)], fill.by = "ident",cols = cols2[!is.na(names(cols2))],stack = T, flip = T,  split.by  = "orig.ident",features =c("Tnf","Il1a","Il1b","Il6","Nlrp3"),slot = "data", assay = "RNA", pt.size = 0.2)+ggtitle("Pro-inflammatory")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Revision/ProInflammatory_VlnPlot.pdf", height = 6, width = 7)


#Immunosup
VlnPlot(RasVCtl_stages.combined,idents =levels(RasVCtl_stages.combined)[c(1,3,2,4:5,7:8)],fill.by = "ident",cols = cols2[!is.na(names(cols2))],stack = T, flip = T,  split.by  = "orig.ident",features =c("Arg1","Il10","Nos2","Ido1","Vegfa","Tgfb1"),slot = "data", assay = "RNA", pt.size = 0.2)+ggtitle("Immunosuppressive")


#Pro Inflammatory + immunosuppressive
VlnPlot(RasVCtl_stages.combined,idents =levels(RasVCtl_stages.combined)[c(1,3,2,4:5,7:8)], fill.by = "ident",cols = cols2[!is.na(names(cols2))],stack = T, flip = T,  split.by  = "orig.ident",features =c("Tnf","Il1a","Il1b","Nlrp3","Arg1","Nos2","Vegfa","Tgfb1"),slot = "data", assay = "RNA", pt.size = 0.2)+ggtitle("Pro-inflammatory")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Revision/ProInflammatory_Immunosup_VlnPlot.pdf", height = 10, width = 7)

#MDSCs FROM LITERATURE:PMN
VlnPlot(RasVCtl_stages.combined,fill.by = "ident",stack = T, flip = T,  features =c("S100a8","S100a9","Anxa1","Lyz2","Csf1","Arg1","Il1b","Cd36","Tgfb1","Cd244a","Vegfa","Stat1","Stat6","Ptgs2","Irf1","Slc27a2"),slot = "data", assay = "RNA", pt.size = 0.2)+ggtitle("PMN-MDSCs")+NoLegend()

#MDSCs FROM LITERATURE: MONO
VlnPlot(RasVCtl_stages.combined,fill.by = "ident",stack = T, flip = T,  features =c("S100a8","S100a9","Arg1","Nos2","Il10","Vegfa","Tnf","Stat3","Tgfb1","Cd274","Tnfrsf10b"),slot = "data", assay = "RNA", pt.size = 0.2)+ggtitle("M-MDSCs")+NoLegend()

#Endothelial
VlnPlot(RasVCtl_stages.combined, features =c("Ptprb","Stab2","Selenop","Igfbp7","Kdr"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Endothelial cells")+NoLegend()

#Lympoid: Granulocytes 
FeaturePlot(RasVCtl_stages.combined, features = c("Adgre1","Ly6g","Ccr3","Siglecf","Fcer1a","Kit","Fcer2a")) 
VlnPlot(RasVCtl_stages.combined, features =c("Adgre1","Ly6g","Ccr3","Siglecf","Fcer1a","Kit","Fcer2a"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Granulocytes")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Cd14_Cd14_umapSCT.pdf")


#Lympoid: NKT 
FeaturePlot(RasVCtl_stages.combined, features = c("Cd3e","Klrb1b","Klrb1c")) 
VlnPlot(RasVCtl_stages.combined, features = c("Cd3e","Klrb1b","Klrb1c"),stack = T,flip = T, pt.size = 0.2)+ggtitle("NKT cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Cd14_Cd14_umapSCT.pdf")

#Myeloid: Neutrophils cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Arg1","Ly6c1","Ly6g")) 
VlnPlot(RasVCtl_stages.combined, features = c("Camp","Ngp","Prtn3","Elane","Ltf"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Neutrophils")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/MDSCcells_VlnplotSCT.pdf")


#Myeloid: gMDSC cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Arg1","Ly6c1","Ly6g")) 
VlnPlot(RasVCtl_stages.combined, features = c("S100a8","S100a9","Cxcl2","Cxcr2","Il1r2"),stack = T,flip = T, pt.size = 0.2)+ggtitle("gMDSC cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/gMDSCcells_VlnplotSCT.pdf")

#Myeloid: moMDSC cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Arg1","Ly6c1","Ly6g")) 
VlnPlot(RasVCtl_stages.combined, features = c("C1qa","Arg1","C1qb","C1qc","Csf1r","S100a4"),stack = T,flip = T, pt.size = 0.2)+ggtitle("moMDSC cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/MDSCcells_VlnplotSCT.pdf")

#Myeloid: pDCs cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Bst2","Siglech","Cd83")) 
VlnPlot(RasVCtl_stages.combined, features = c("Itgam","Bst2","Siglech","Cd83"),stack = T,flip = T, pt.size = 0.2)+ggtitle("pDCs cells")
VlnPlot(RasVCtl_stages.combined, features = c("Grm8","Tcf4","Bst2","Runx2","Irf8","Siglech"),stack = T,flip = T, pt.size = 0.2)+ggtitle("pDCs cells")

ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/pDCscells_VlnplotSCT.pdf")


#Myeloid:cDCs cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Itgax","H2-Ab1","Xcr1","Clec9a","Cd8a","Itgae","Sirpa","Cd83")) 
VlnPlot(RasVCtl_stages.combined, features = c("Itgam","Itgax","H2-Ab1","Xcr1","Clec9a","Cd8a","Itgae","Sirpa","Cd83"),stack = T,flip = T, pt.size = 0.2)+ggtitle("cDCs cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/cDCscells_VlnplotSCT.pdf")

#Myeloid:cDC1s cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Itgax","H2-Ab1","Xcr1","Clec9a","Cd8a","Itgae","Sirpa","Cd83")) 
VlnPlot(RasVCtl_stages.combined, features = c("Cst3","Ppt1","H2-Ab1","Xcr1","Clec9a","Flt3"),stack = T,flip = T, pt.size = 0.2)+ggtitle("cDC1s cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/cDCscells_VlnplotSCT.pdf")

#Myeloid:Macrophages cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Adgre1","Cd86","Cd80","Nos2","Mrc1","Arg1","Clec4f","Tlr4","Timd4")) 
VlnPlot(RasVCtl_stages.combined, features = c("Spi1","Itgam","Adgre1","Cd86","Cd80","Nos2","Mrc1","Arg1","Clec4f","Tlr4","Timd4","Vsig4"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Macrophages cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Macrophages_VlnplotSCT.pdf")

#Myeloid:Monocytes cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Itgam","Cd14","Ccr2","Cx3cr1","Ly6c1")) 
VlnPlot(RasVCtl_stages.combined, features = c("Itgam","Cd14","Ccr2","Cx3cr1","Ly6c1"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Monocytes")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Monocytes_VlnplotSCT.pdf")


#Lympoid: T cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Cd3e","Tcrac","Tcrb","Cd8a","Cd4","Il2ra","Cd69","Gzmb","Prf1","Cd44","Il17ra","Pdcd1","Foxp3")) 
VlnPlot(RasVCtl_stages.combined, features = c("Trbc2","Cd3d","Cd3g","Trbc1","Cd3e","Tcf7"),stack = T,flip = T, pt.size = 0.2)+ggtitle("T cells")
VlnPlot(RasVCtl_stages.combined, features = c("Cd3e","Trac","Cd8a","Cd4","Il2ra","Cd69","Gzmb","Prf1","Cd44","Il17ra","Pdcd1","Ctla4","Foxp3"),stack = T,flip = T, pt.size = 0.2)+ggtitle("T cells")

ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Tcells_VlnplotSCT.pdf")

#Lympoid: NK cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Klrb1b","Klrb1c","Ncr1","Klrk1","Cd69","Gzmb","Prf1")) 
VlnPlot(RasVCtl_stages.combined, features = c("Klrb1b","Klrb1c","Ncr1","Klrk1","Cd69","Gzmb","Prf1"),stack = T,flip = T, pt.size = 0.2)+ggtitle("NK cells")
VlnPlot(RasVCtl_stages.combined, features = c("Gzma","Xcl1","Gzmb","Klrb1c","Klrd1","Klrk1"),stack = T,flip = T, pt.size = 0.2)+ggtitle("NK cells")

ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/NKcells_VlnplotSCT.pdf")

FeaturePlot(RasVCtl_stages.combined, features = c("Lgals3","Plin2","Abca1","Trem2","Nr1h3","Irf4","Dhcr24","Gpnmb","Ccl6")) 

#Lympoid: NKT cells 
VlnPlot(RasVCtl_stages.combined, features = c("St6galnac3","Cd226","Kcnq5","Skap1","Atp8a2","Inpp4b"),stack = T,flip = T, pt.size = 0.2)+ggtitle("NK cells")

#Lympoid: B cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Cd19","Ighd","Cd27")) 
VlnPlot(RasVCtl_stages.combined, features = c("Cd19","Ighd","Cd27","Ighm"),stack = T,flip = T, pt.size = 0.2)+ggtitle("B cells")
VlnPlot(RasVCtl_stages.combined, features = c("Igkc","Bank1","Ighm","Cd79b","Cd79a"),stack = T,flip = T, pt.size = 0.2)+ggtitle("B cells")

ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Bcells_VlnplotSCT.pdf")

#Lympoid: Plasma cells 
FeaturePlot(RasVCtl_stages.combined, features = c("Tnfrsf17","Sdc1")) 
VlnPlot(RasVCtl_stages.combined, features = c("Tnfrsf17","Sdc1","Mzb1"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Plasma cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Plasmacells_VlnplotSCT.pdf")


#CD45neg: 
c("Fxyd2","Tff3","Sstr2",#cholangiocytes
  "Col6a2","Pdgfrb","Vim", #Fibroblasts
  "Pecam1","Egfl7",#Endothelial
  "Col1a1","Igfbp3","Clec4g",#HSC
  "C1qa","Cd163",#KC
  "Alb","Cyp2e1")#hepa
VlnPlot(RasVCtl_stages.combined, features = c("Fxyd2","Tff3","Sstr2",
                                                      "Col6a2","Pdgfrb",
                                                      "Pecam1","Egfl7",
                                                      "Col1a1","Igfbp3","Col1a2","Clec4g",
                                                      "C1qa","Cd163",
                                                      "Alb","Cyp2e1"),stack = T,flip = T, pt.size = 0.2)+ggtitle("Plasma cells")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/CD45neg_VlnplotSCT.pdf")


#Heatmap of cell type markers
library(dplyr)

RasVCtl_stages.combined.markers[RasVCtl_stages.combined.markers$gene %in%c("S100a8","S100a9","Il1b", "Arg1","Lyz2","Csf1r","Igkc" ,"Cd79a","Ighm","Trbc2","Cd3d" ,"Tcf7","Gzma","Xcl1" ,"Klrb1c","Cd226","St6galnac3","Kcnq5", "Bst2","Ly6d" ,"Grm8","Camp" , "Prtn3","Elane" ,"Xcr1","Cst3","H2-Ab1"),]%>%
  group_by(cluster) %>% top_n(5, avg_log2FC) -> top
DoHeatmap(object =  RasVCtl_stages.combinedminusEndothelial,features = top$gene, label = T,group.colors=cols[!is.na(names(cols))])
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/CelltypeDEGmarkerspercluster_5genes.pdf", height = 5, width = 12)

############################## 


############################### 9. Differential expression ############################### 
############################## 
# For performing differential expression after integration, we switch back to the original data
# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(RasVCtl_stages.combined) <- "RNA"

library(future)

# change the current plan to access parallelization
plan("multisession", workers = 6)

RasVCtl_stages.combined.markers <- FindAllMarkers(RasVCtl_stages.combined, logfc.threshold = 0.5,  min.pct =  0.25, 
                                                          min.diff.pct = 0.25,assay = "RNA",slot="data", features = features)

future:::ClusterRegistry("stop")

t<-RasVCtl_stages.combined.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]#0]
c1<-t$gene[t$cluster==1]#1]
c2<-t$gene[t$cluster==2]#2]
c3<-t$gene[t$cluster==3]#3]
c4<-t$gene[t$cluster==4]#4]
c5<-t$gene[t$cluster==5]#5]
c6<-t$gene[t$cluster==6]#6]
c7<-t$gene[t$cluster==7]#7]
c8<-t$gene[t$cluster==8]
c9<-t$gene[t$cluster==9]
c10<-t$gene[t$cluster==10]
c11<-t$gene[t$cluster==11]
c12<-t$gene[t$cluster==12]
#c13<-t$gene[t$cluster==13]
#c14<-t$gene[t$cluster==14]
#c15<-t$gene[t$cluster==15]
#c16<-t$gene[t$cluster==16]
#c17<-t$gene[t$cluster==17]
#c18<-t$gene[t$cluster==18]
#c19<-t$gene[t$cluster==19]
markerseachcluster<-data.frame(c0,c1, c2, c3, c4, c5, c6, c7,c8,c9,c10)
###########Export as excel
library(openxlsx)
SourceData <- createWorkbook()

addWorksheet(SourceData, "markerseachcluster_top50_RasV")#Add sheets

# Write the data to the sheets
writeData(SourceData, sheet = "markerseachcluster_top50_RasV", x = markerseachcluster )

# Export the file
saveWorkbook(SourceData, "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/DEG_MarkersClusters_50top_RasV_res025.xlsx")

RasVCtl_stages.combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

#DefaultAssay(RasVCtl_stages.combined) <- "integrated"
#DoHeatmap(object =  RasVCtl_stages.combined,features = top10$gene, label = T,group.colors=cols[!is.na(names(cols))]) + NoLegend()#+ scale_fill_gradientn(colors = ArchRPalettes$coolwarm)

#RasVCtl_stages.combined.markers[RasVCtl_stages.combined.markers$gene %in%c("Tnf","Il1a","Il1b","Il6","Nlrp3","Arg1","Il10","Nos2","Ido1","Vegfc","Tgfb1","Vegfa","Vegf","Tgfb2","Cd244a","Vegfa","Ptgs2","Arg1","Nos2","Il10","Vegfa","Tnf","Stat3","Tgfb1","Cd274","Tnfrsf10b"),]%>%
#  group_by(cluster) -> top
DoHeatmap(object =  RasVCtl_stages.combined,features = top10$gene, label = T,group.colors=cols[!is.na(names(cols))])
#ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/CelltypeDEGimmunesuppinflgenes.pdf", height = 5, width = 12)
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Revision/Heatmap_CelltypeDEGtop10.pdf", height = 5, width = 12)

######Geneset enrichment
#Try enrichment 
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)

#Keep all the DEG
t<-RasVCtl_stages.combined.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05) 
c0<-t$gene[t$avg_log2FC>0 &t$cluster=="Monocytic cells"]#0]
c1<-t$gene[t$avg_log2FC>0 &t$cluster=="Neutrophils"]#1]
c2<-t$gene[t$avg_log2FC>0 &t$cluster=="T cells"]#2]
c3<-t$gene[t$avg_log2FC>0 &t$cluster=="cDC1"]#3]
c4<-t$gene[t$avg_log2FC>0 &t$cluster=="B cells"]#4]
c5<-t$gene[t$avg_log2FC>0 &t$cluster=="migratory DCs"]#5]
c6<-t$gene[t$avg_log2FC>0 &t$cluster=="Platelets"]#6]
c7<-t$gene[t$avg_log2FC>0 &t$cluster=="Endothelial cells"]#7]
c8<-t$gene[t$avg_log2FC>0 &t$cluster=="pDC"]
#c9<-t$gene[t$avg_log2FC>0 &t$cluster==9]
#RasVCtl_stages.combined.markers_topUPDEGList<-list(c0,c1,  c2,  c3,c4,c5,c6,c7,c8,c9)
RasVCtl_stages.combined.markers_topUPDEGList<-list(c0,c1,  c2,  c3,c4,c5,c6,c7,c8)

#Conver genename to entrez id
geneid.ls <- RasVCtl_stages.combined.markers_topUPDEGList %>% map(~{
  
  
  gene.df <- select(org.Mm.eg.db,
                    keys = .x,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})

kegg.ls <- geneid.ls %>% map(~{
  eKEGG <- enrichKEGG(
    gene = .x,
    pvalueCutoff = 0.01, 
    organism = 'mmu',
    pAdjustMethod = "BH"
  )
  eKEGG<-setReadable(eKEGG, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  return(eKEGG)
})
names(geneid.ls)<-levels(Idents(RasVCtl_stages.combined))
ck <- compareCluster(geneCluster = geneid.ls, fun = enrichKEGG,pvalueCutoff = 0.01, 
                     #OrgDb = org.Mm.eg.db,
                     organism = 'mmu',
                     
                     pAdjustMethod = "BH")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 
tail(ck)
dotplot(ck)+theme_bw()+
  theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.justification=c(1,0), 
    #legend.title = element_text("Clusters"),  
    axis.text.x  = element_text(angle=45, hjust = 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Enrichment_RasDandCTLallclusters_KEGG.pdf", width = 9, height = 12)




############################### 


############################### 11. Assign cell annot ############################### 
############################### 

new.cluster.ids <- c("Endothelial cells", "T cells","B cells", "Monocytic cells",  "NK cells", 
                     "Neutrophils", "Hepatocytes", "pDCs","Plasma cells","Hepatocytes II","pDCs","Fibroblasts")
names(new.cluster.ids) <- levels(RasVCtl_stages.combined)
RasVCtl_stages.combined <- RenameIdents(RasVCtl_stages.combined, new.cluster.ids)
DimPlot(RasVCtl_stages.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

RasVCtl_stages.combined$Idents<-Idents(RasVCtl_stages.combined)
############################### 


############################### 12. Proportion of cell types per condition ############################### 
############################### 
RasVCtl_stages.combined$idents<-Idents(RasVCtl_stages.combined)
proportions<-RasVCtl_stages.combined@meta.data %>% group_by(sample_name,idents, replicate) %>% summarise(sample_name,idents, replicate) %>% table()

#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:3, function(x){proportions[x,,replicate = "R1"]/t[x]})

colnames(t2)<-rownames(proportions)
#rownames(t2)<-new.cluster.ids
#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:3, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)

#R3
t5<-apply(proportions[,,replicate = "R3"], 1, sum)
t6<-sapply(1:3, function(x){proportions[x,,replicate = "R3"]/t5[x]})
colnames(t6)<-rownames(proportions)

t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)),
         data.frame("rep"="R3",reshape2::melt(t6)))
colnames(t)<-c("Rep","cluster","Model","value")
ggplot(t, aes(x=Model,y=value, fill=(cluster)) )+geom_bar(stat="identity", col="white")+theme_bw()+
  scale_fill_manual(values =cols[!is.na(names(cols))])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        #legend.title = element_text("Clusters"),  
        axis.text.x  = element_text(angle=0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))


names(cols)<-names(table(Idents(RasVCtl_stages.combined)))
t$Model
library(ggpubr)
library(rstatix)

t$Model<-factor(t$Model, levels = c("Ctl","RasV_Pretumor","RasV_Endpoint"))

#t$type<-NA
#t$type[t$cluster%in% c("Neutrophils","Monocytic cells","cDC1","pDCs","migratory DCs")]<-"Myeloid"
#t$type[t$cluster%in% c("T cells","Naive T cells","B cells","NK cells","Plasma cells")]<-"Lymphoid"
#t$type[t$cluster%in% c("Endothelial cells","HSC+Hepatocytes","Mixed(Hepatocytes+Fibro)","Hepatocytes")]<-"Other"

p <- ggbarplot(t[complete.cases(t),], x = "Model", y = "value", add = "mean_se",scales = "free_y",
               
               fill =  "cluster",palette = cols[!is.na(names(cols))])+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x  = element_text(angle=45, hjust=1))
p

#
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Proportions_annotatedclusted_dividebytype_forB3B6.pdf", width = 4, height = 6)

#Comparison between controls:
proportions<-RasVCtl_stages.combined@meta.data %>% filter(orig.ident=="Ctl") %>% group_by(stage, idents, replicate) %>% summarise(stage,idents, replicate) %>% table()

#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:2, function(x){proportions[x,,replicate = "R1"]/t[x]})

colnames(t2)<-rownames(proportions)
#rownames(t2)<-new.cluster.ids
#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:2, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)


t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)))
colnames(t)<-c("Rep","cluster","Model","value")


names(cols)<-names(table(Idents(RasVCtl_stages.combined)))
t$Model
library(ggpubr)
library(rstatix)


p <- ggbarplot(t[t$Model%in%c("Intermediate","OLD_Intermediate"),], x = "Model", y = "value", add = "mean_se",scales = "free_y",
               
               fill =  "cluster",palette = cols[!is.na(names(cols))])+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x  = element_text(angle=45, hjust=1))
p

#
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Proportions_annotatedclusted_dividebytype.pdf", width = 8, height = 6)

############################### 


############################### 13. Differential abundance:MiloR############################### 
############################### 
library(miloR)
library(SingleCellExperiment)

#1. Create Milo object
DefaultAssay(RasVCtl_stages.combined) <- "integrated"
t<-subset(RasVCtl_stages.combined, subset=sample_name %in% c("Ctl","RasV_Pretumor"))
RasVCtl_stages.combined_sce <- as.SingleCellExperiment(t)
RasVCtl_stages.combined_milo <- Milo(RasVCtl_stages.combined_sce)

#2. Visualize the data
plotReducedDim(RasVCtl_stages.combined_milo, colour_by="ident", dimred = "UMAP") 

#3. Differential abundance testing
#3.2 Construct KNN graph
traj_milo<-buildGraph(RasVCtl_stages.combined_milo, k = 30, d = 50)

#3.3 Defining representative neighbourhoods on the KNN graph
traj_milo <- makeNhoods(traj_milo, prop = 0.075, k = 30, d=50, refined = TRUE)
plotNhoodSizeHist(traj_milo)

#3.4 Counting cells in neighbourhoods
meta.df <- t@meta.data
#meta.df$SampID <- paste(meta.df$orig.ident, Idents(RasVCtl_stages.combined), sep="_")
#meta.df$Ident<-Idents(RasVCtl_stages.combined)
#traj_milo$SampID <- paste(traj_milo$orig.ident, meta.df$Ident, sep="_")
meta.df$replicate[meta.df$stage=="OLD_Intermediate"]<-"R3"#change the number of the replicate for the OLD control
meta.df$SampID <- paste(meta.df$sample_name, meta.df$replicate, sep="_")
meta.df$SampID2 <- paste(meta.df$SampID, meta.df$idents, sep="_")
#meta.df$Ident<-Idents(RasVCtl_stages.combined)
traj_milo$SampID <- meta.df$SampID




traj_milo <- countCells(traj_milo, meta.data=meta.df, samples="SampID")
head(nhoodCounts(traj_milo))

#3.5 Defining experimental design
test.meta <- data.frame("SampID"=names(table(meta.df$SampID)), "orig.ident"=c(rep("Ctl", 3),rep("RasV_Pretumor", 3)))

#test.meta <- data.frame("SampID"=names(table(meta.df$SampID)), "orig.ident"=c(rep("GMCSF", 9),rep("IgG", 9)))
#test.meta$Sample <- paste(test.meta$orig.ident, test.meta$Replicate, sep="_")
rownames(test.meta) <- test.meta$SampID

#Doesnt work, need replicates
#test.meta <- data.frame("orig.ident"=names(table(meta.df$orig.ident)), "replicate"=names(table(meta.df$replicate)))
#rownames(test.meta) <- test.meta$orig.ident


#3.6 Computing neighbourhood connectivity
traj_milo <- calcNhoodDistance(traj_milo, d=50)

#3.7 Testing
da.res <- testNhoods(traj_milo, design=~orig.ident, design.df=test.meta[colnames(nhoodCounts(traj_milo)), ])

da.res %>%
  arrange(SpatialFDR) %>%
  head() 

#3.8 Inspecting DA testing results
ggplot(da.res, aes(PValue)) + geom_histogram(bins=50)
ggplot(da.res, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) 


traj_milo <- buildNhoodGraph(traj_milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(traj_milo, dimred = "UMAP", colour_by="ident", text_by = "ident", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da.res, layout="UMAP",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_ControlvsPretumor_MiloR.pdf", width = 11, height = 5)

da.res <- annotateNhoods(traj_milo, da.res, coldata_col = "ident")
head(da.res)

ggplot(da.res, aes(ident_fraction)) + geom_histogram(bins=50)

da.res$ident <- ifelse(da.res$ident_fraction < 0.7, "Mixed", da.res$ident)

#da.res$ident<-factor(da.res$ident, levels=rev(c(levels(RasVCtl_stages.combined), "Mixed")))

plotDAbeeswarm(da.res, group.by = "ident")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/plotDAbeeswarm_Ras12VendpointvsPretumor_MiloR.pdf", width = 5, height = 6)



############################### 

############################### 14. Zooim into cluster############################### 
#-----------------Cluster 0 : Endothelial cells-----------------
clust0_Endothelialcells<-subset(x = RasVCtl_stages.combined, idents = c("Endothelial cells"))
clust0_Endothelialcells

DefaultAssay(clust0_Endothelialcells) <- "integrated"
clust0_Endothelialcells <- RunPCA(clust0_Endothelialcells, verbose = T)
DimPlot(clust0_Endothelialcells, reduction = "pca")
ElbowPlot(clust0_Endothelialcells, reduction = "pca",ndims = 50)

clust0_Endothelialcells <- RunUMAP(clust0_Endothelialcells, reduction = "pca", dims = 1:50)
clust0_Endothelialcells <- FindNeighbors(clust0_Endothelialcells, reduction = "pca", dims = 1:50)
clust0_Endothelialcells <- FindClusters(clust0_Endothelialcells, resolution = 0.07)

#Assign new identities: based on both DEG and transfer guilliams
new.cluster.ids <- c("moDCs", "Classical monocytes", "BMDMs", "cDC2s","Proliferating cDC1s","Inflammatory Monocytes", "Patrolling Monocytes", "KCs",
                     "Migratory DCs")
names(new.cluster.ids) <- levels(clust0_Endothelialcells)
clust0_Endothelialcells <- RenameIdents(clust0_Endothelialcells, new.cluster.ids)

DefaultAssay(clust0_Endothelialcells) <- "SCT"
FeaturePlot(clust0_Endothelialcells,features =  c("Flt1"))

# Visualization
# Visualization
#cols_clust0_Endothelialcells<-sample(x = ArchRPalettes[[3]], size = 8, replace = F)
names(cols_clust0_Endothelialcells)<-names(table(Idents(clust0_Endothelialcells)))

p1 <- DimPlot(clust0_Endothelialcells, reduction = "umap", cols = cols2[!is.na(names(cols2))], group.by = "orig.ident")#, cells = cDC1scells)+ylim(c(-15,8))+xlim(c(-5,10))#,group.by = "orig.ident")
p2 <- DimPlot(clust0_Endothelialcells, reduction = "umap",cols = cols_clust0_Endothelialcells[!is.na(names(cols_clust0_Endothelialcells))],label = TRUE, repel = TRUE)
p1+p2
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_controlRas12V_zoomC0_Endotheliacells_res05_annotatedclusted.pdf", width = 10, height = 5)

save(clust0_Endothelialcells, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust0_Endothelialcells.RData")
#load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust0_Endothelialcells.RData")


#DEG
DefaultAssay(clust0_Endothelialcells) <- "RNA"

plan("multisession", workers = 6)
clust0_Endothelialcells.markers <- FindAllMarkers(clust0_Endothelialcells, logfc.threshold = 0.5,  min.pct =  0.25, 
                                                  min.diff.pct = 0.25,assay = "RNA",slot="data")
future:::ClusterRegistry("stop")


t<-clust0_Endothelialcells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]
c5<-t$gene[t$cluster==5]
c6<-t$gene[t$cluster==6]
c7<-t$gene[t$cluster==7]
c8<-t$gene[t$cluster==8]

clust0_Endothelialcells.markers_top50DEG<-cbindX( data.frame(c0), data.frame(c1),  data.frame(c2),  data.frame(c3),  data.frame(c4))#,  data.frame(c5),  data.frame(c6),  data.frame(c7))


###########Export as excel
library(openxlsx)
SourceData <- createWorkbook()

addWorksheet(SourceData, "markerseachcluster_c0_Endo")#Add sheets

# Write the data to the sheets
writeData(SourceData, sheet = "markerseachcluster_c0_Endo", x = clust0_Endothelialcells.markers_top50DEG )

# Export the file
saveWorkbook(SourceData, "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominC0_Endothelialcells_DEG_MarkersClusters_RasV_res05.xlsx")






#Proportions
library(plyr)
proportions<-clust0_Endothelialcells@meta.data %>% group_by(sample_name,predicted.id, replicate) %>% summarise(orig.ident,predicted.id, replicate) %>% table()
proportions<-clust0_Endothelialcells@meta.data %>% group_by(sample_name,Idents(clust0_Endothelialcells), replicate) %>% plyr:::summarise(sample_name,Idents(clust0_Endothelialcells), replicate) %>% table()

#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:3, function(x){proportions[x,,replicate = "R1"]/t[x]})
colnames(t2)<-rownames(proportions)
#rownames(t2)<-new.cluster.ids

#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:3, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)

#R3
t5<-apply(proportions[,,replicate = "R3"], 1, sum)
t6<-sapply(1:3, function(x){proportions[x,,replicate = "R3"]/t5[x]})
colnames(t6)<-rownames(proportions)

t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)),
         data.frame("rep"="R3",reshape2::melt(t6)))
colnames(t)<-c("Rep","cluster","Model","value")


ggplot(t, aes(x=Model,y=value, fill=as.factor(cluster)) )+geom_bar(stat="identity", col="white")+theme_bw()+
  scale_fill_manual(values=cols4[!is.na(names(cols4))])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        #legend.title = element_text("Clusters"),  
        axis.text.x  = element_text(angle=0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Proportions_controlRas12D_zoomcluster19_cluster.pdf", width = 3.5, height = 6)


library(ggpubr)
library(rstatix)

t$Model<-factor(t$Model, levels = c("Ctl","RasV_Pretumor","RasV_Endpoint"))

t$type<-NA
t$type[t$cluster%in% c("Neutrophils","Monocytic cells","cDC1","pDCs","migratory DCs")]<-"Myeloid"
t$type[t$cluster%in% c("T cells","Naive T cells","B cells","NK cells","Plasma cells")]<-"Lymphoid"
t$type[t$cluster%in% c("Endothelial cells","HSC+Hepatocytes","Mixed(Hepatocytes+Fibro)","Hepatocytes")]<-"Other"

t$cluster<-as.character(t$cluster)
p <- ggbarplot(t[complete.cases(t),], x = "Model", y = "value", add = "mean_se",
               fill =  "cluster",palette = cols_clust0_Endothelialcells[!is.na(names(cols_clust0_Endothelialcells))])+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x  = element_text(angle=45, hjust=1))
p

#
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Proportions_ZoominC0_Endothelialcells.pdf", width = 4, height = 6)


#Try enrichment 
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)

#Keep all the DEG
t<-clust14_TcellsNaive.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05) 
c0<-t$gene[t$avg_log2FC>0 &t$cluster==0]
c1<-t$gene[t$avg_log2FC>0 &t$cluster==1]
c2<-t$gene[t$avg_log2FC>0 &t$cluster==2]
c3<-t$gene[t$avg_log2FC>0 &t$cluster==3]
c4<-t$gene[t$avg_log2FC>0 &t$cluster==4]
c5<-t$gene[t$avg_log2FC>0 &t$cluster==5]
c6<-t$gene[t$avg_log2FC>0 &t$cluster==6]
c7<-t$gene[t$avg_log2FC>0 &t$cluster==7]
c8<-t$gene[t$avg_log2FC>0 &t$cluster==8]

clust14_TcellsNaive.markers_topUPDEGList<-list(c0,c1,  c2,  c3,c4,c5)

#Conver genename to entrez id
geneid.ls <- clust14_TcellsNaive.markers_topUPDEGList %>% map(~{
  
  
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})


names(geneid.ls)<-levels(Idents(clust14_TcellsNaive))

#Try hallmarks msigdb: 
library(msigdbr)
m_df <- msigdbr(species = "Mus musculus")
table(m_df$gs_subcat)
m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
ck <- compareCluster(geneCluster = geneid.ls, fun = enricher,TERM2GENE=m_t2g,pvalueCutoff = 0.05, 
                     #OrgDb = org.Mm.eg.db,
                     #organism = 'mmu',
                     
                     pAdjustMethod = "BH")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 
tail(ck)
#showCategory=NULL
dotplot(ck)+theme_bw()+
  theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.justification=c(1,0), 
    #legend.title = element_text("Clusters"),  
    axis.text.x  = element_text(angle=45, hjust = 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),text=element_text(family="Helvetica", size=8, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Enrichment_C19_GOBP.pdf", width = 6, height = 10)











#-----------------Cluster 1 : T cells-----------------
clust1_Tcells<-subset(x = RasVCtl_stages.combined, idents = c("T cells"))
clust1_Tcells

DefaultAssay(clust1_Tcells) <- "integrated"
clust1_Tcells <- RunPCA(clust1_Tcells, verbose = T)
DimPlot(clust1_Tcells, reduction = "pca")
ElbowPlot(clust1_Tcells, reduction = "pca",ndims = 50)

clust1_Tcells <- RunUMAP(clust1_Tcells, reduction = "pca", dims = 1:50)
clust1_Tcells <- FindNeighbors(clust1_Tcells, reduction = "pca", dims = 1:50)
clust1_Tcells <- FindClusters(clust1_Tcells, resolution = 0.05)

#Assign new identities: based on both DEG and transfer guilliams
new.cluster.ids <- c("moDCs", "Classical monocytes", "BMDMs", "cDC2s","Proliferating cDC1s","Inflammatory Monocytes", "Patrolling Monocytes", "KCs",
                     "Migratory DCs")
names(new.cluster.ids) <- levels(clust1_Tcells)
clust1_Tcells <- RenameIdents(clust1_Tcells, new.cluster.ids)

new.cluster.ids <- c("C1_Monocytic_March3+", "C2_Classical monocytes_Chil3+", "C3_Monocytic_C1qa+", "C4_cDC2_Cd209a+","C5_cDC1s_Wdfy4+", "C6_Monocytic_Cxcl3+", "C7_Monocytic_Adgre4+",
                     "C8_KC_Slc40a1+", "C9_Migratory DCs_Ccr7+")
names(new.cluster.ids) <- levels(clust1_Tcells)
clust1_Tcells <- RenameIdents(clust1_Tcells, new.cluster.ids)

# Visualization
# Visualization
#cols4<-sample(ArchRPalettes[[6]], size = 8, replace = F)
names(cols4)<-names(table(Idents(clust1_Tcells)))
cols<-ArchRPalettes[[8]][c(3,10,11)]
names(cols)<-names(table(clust1_Tcells$sample_name))

p1 <- DimPlot(clust1_Tcells, reduction = "umap", cols = cols[!is.na(names(cols))], group.by = "sample_name", shuffle = T)#, cells = cDC1scells)+ylim(c(-15,8))+xlim(c(-5,10))#,group.by = "orig.ident")
p2 <- DimPlot(clust1_Tcells, reduction = "umap",cols = cols4[!is.na(names(cols4))],label = TRUE, repel = TRUE, shuffle = T)
p1+p2
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_controlRas12V_zoomC1_Tcells_res05_unannotatedclusted_forB3B6.pdf", width = 10, height = 5)

save(clust1_Tcells, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust1_Tcells.RData")
load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust1_Tcells.RData")

FeaturePlot(clust1_Tcells,features =  c("Pdcd1","Cd274","Ctla4"))
#DEG
DefaultAssay(clust1_Tcells) <- "RNA"

plan("multisession", workers = 6)
clust1_Tcells.markers <- FindAllMarkers(clust1_Tcells, logfc.threshold = 0.5,  min.pct =  0.25, 
                                               min.diff.pct = 0.25,assay = "RNA",slot="data")
future:::ClusterRegistry("stop")


t<-clust1_Tcells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]
c5<-t$gene[t$cluster==5]
c6<-t$gene[t$cluster==6]
c7<-t$gene[t$cluster==7]
c8<-t$gene[t$cluster==8]
c9<-t$gene[t$cluster==9]
c10<-t$gene[t$cluster==10]
c11<-t$gene[t$cluster==11]

clust1_Tcells.markers_top50DEG<-cbindX( data.frame(c0), data.frame(c1),  data.frame(c2),  data.frame(c3),  data.frame(c4),  data.frame(c5),  data.frame(c6),  data.frame(c7))#,  data.frame(c8),  data.frame(c9),  data.frame(c10),  data.frame(c11))

DefaultAssay(clust1_Tcells) <- "SCT"
FeaturePlot(clust1_Tcells, features = c("Csf2ra","Csf2rb","Ccl6","Ccl17","Ccl24")) #C
VlnPlot(clust1_Tcells, features = c("Il23a"), group.by = "orig.ident") #C
library(readxl)
Tcellsignatures_public<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/Public/Tcellsignatures.xlsx")

Tcellsignatures_public_Exhaustion<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion[!is.na(Tcellsignatures_public$Exhaustion)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_G1S<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G1/S`[!is.na(Tcellsignatures_public$`G1/S`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_G2M<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G2/M`[!is.na(Tcellsignatures_public$`G2/M`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_Exhaustion_CD8<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD8[!is.na(Tcellsignatures_public$Exhaustion_CD8)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_Exhaustion_CD4<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD4[!is.na(Tcellsignatures_public$Exhaustion_CD4)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)


Tcellsignatures_public_list<-list("Memory"=Tcellsignatures_public$Memory,"Effector"=Tcellsignatures_public$Effector,
                                  "Exhaustion"=Tcellsignatures_public_Exhaustion$MGI.symbol,
                                  "G1S"=Tcellsignatures_public_G1S$MGI.symbol,  "G2M"=Tcellsignatures_public_G2M$MGI.symbol,
                                  "Exhaustion_CD8"=Tcellsignatures_public_Exhaustion_CD8$MGI.symbol,"Exhaustion_CD4"=Tcellsignatures_public_Exhaustion_CD4$MGI.symbol,
                                  "Progenitor_Exh"=Tcellsignatures_public$Progenitor_Exh,"Effector_like"=Tcellsignatures_public$Effector_like,
                                  "Terminally_Exh"=Tcellsignatures_public$Terminally_Exh, "Proliferating"=Tcellsignatures_public$Proliferating)
clust1_Tcells <- AddModuleScore(
  object = clust1_Tcells,
  features = Tcellsignatures_public_list,
  ctrl = 5,
  name = 'Tcellsignatures'
)
library(RColorBrewer)
p1<-FeaturePlot(object = clust1_Tcells,  min.cutoff = -2, max.cutoff = 5,features = c("Tcellsignatures1"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Memory")
p2<-FeaturePlot(object = clust1_Tcells, min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures2"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Effector")
p3<-FeaturePlot(object = clust1_Tcells,min.cutoff = -5, max.cutoff = 5, features = c("Tcellsignatures3"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Exhaustion")
p4<-FeaturePlot(object = clust1_Tcells,min.cutoff = -1, max.cutoff = 1, features = c("Tcellsignatures4"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("G1S")
p5<-FeaturePlot(object = clust1_Tcells,min.cutoff = -1, max.cutoff = 1, features = c("Tcellsignatures5"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("G2M")
p6<-FeaturePlot(object = clust1_Tcells,  min.cutoff = -3, max.cutoff = 3,features = c("Tcellsignatures6"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Exhaustion_CD8")
p7<-FeaturePlot(object = clust1_Tcells, min.cutoff = -5, max.cutoff = 5, features = c("Tcellsignatures7"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Exhaustion_CD4")
p8<-FeaturePlot(object = clust1_Tcells,min.cutoff = -5, max.cutoff = 5, features = c("Tcellsignatures8"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Progenitor_Exh")
p9<-FeaturePlot(object = clust1_Tcells, min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures9"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Effector_like")
p10<-FeaturePlot(object = clust1_Tcells,min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures10"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Terminally_Exh")
p11<-FeaturePlot(object = clust1_Tcells, min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures11"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Proliferating")

pdf("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominTcells_Tcellsignatures_Enrichment.pdf", width = 13, height = 8)
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, ncol=4)
dev.off()

###########Export as excel
library(openxlsx)
SourceData <- createWorkbook()

addWorksheet(SourceData, "markerseachcluster_c1_Tcells")#Add sheets

# Write the data to the sheets
writeData(SourceData, sheet = "markerseachcluster_c1_Tcells", x = clust1_Tcells.markers_top50DEG )

# Export the file
saveWorkbook(SourceData, "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominC1_Tcells_DEG_MarkersClusters_RasV_res05.xlsx")






#Proportions
library(plyr)
proportions<-clust1_Tcells@meta.data %>% group_by(sample_name,predicted.id, replicate) %>% summarise(orig.ident,predicted.id, replicate) %>% table()
proportions<-clust1_Tcells@meta.data %>% group_by(sample_name,Idents(clust1_Tcells), replicate) %>% plyr:::summarise(sample_name,Idents(clust1_Tcells), replicate) %>% table()

#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:3, function(x){proportions[x,,replicate = "R1"]/t[x]})
colnames(t2)<-rownames(proportions)
#rownames(t2)<-new.cluster.ids

#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:3, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)

#R3
t5<-apply(proportions[,,replicate = "R3"], 1, sum)
t6<-sapply(1:3, function(x){proportions[x,,replicate = "R3"]/t5[x]})
colnames(t6)<-rownames(proportions)

t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)),
         data.frame("rep"="R3",reshape2::melt(t6)))
colnames(t)<-c("Rep","cluster","Model","value")


ggplot(t, aes(x=Model,y=value, fill=as.factor(cluster)) )+geom_bar(stat="identity", col="white")+theme_bw()+
  scale_fill_manual(values=cols4[!is.na(names(cols4))])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        #legend.title = element_text("Clusters"),  
        axis.text.x  = element_text(angle=0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Proportions_controlRas12D_zoomcluster19_cluster.pdf", width = 3.5, height = 6)


library(ggpubr)
library(rstatix)

t$Model<-factor(t$Model, levels = c("Ctl","RasV_Pretumor","RasV_Endpoint"))


t$cluster<-as.character(t$cluster)
p <- ggbarplot(t[complete.cases(t),], x = "Model", y = "value", add = "mean_se",
               
               fill =  "cluster",palette = cols4[!is.na(names(cols4))])+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x  = element_text(angle=45, hjust=1))
p

#
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Proportions_ZoominC1_Tcells_res05.pdf", width = 4, height = 6)


#Try enrichment 
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)

#Keep all the DEG
t<-clust1_Tcells.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05) 
c0<-t$gene[t$avg_log2FC>0 &t$cluster==0]
c1<-t$gene[t$avg_log2FC>0 &t$cluster==1]
c2<-t$gene[t$avg_log2FC>0 &t$cluster==2]
c3<-t$gene[t$avg_log2FC>0 &t$cluster==3]
c4<-t$gene[t$avg_log2FC>0 &t$cluster==4]
c5<-t$gene[t$avg_log2FC>0 &t$cluster==5]
c6<-t$gene[t$avg_log2FC>0 &t$cluster==6]
c7<-t$gene[t$avg_log2FC>0 &t$cluster==7]
c8<-t$gene[t$avg_log2FC>0 &t$cluster==8]

clust1_Tcells.markers_topUPDEGList<-list(c0,c1,  c2,  c3,c4,c5,c6,c7)

#Conver genename to entrez id
geneid.ls <- clust1_Tcells.markers_topUPDEGList %>% map(~{
  
  
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})


names(geneid.ls)<-levels(Idents(clust1_Tcells))

#Try hallmarks msigdb: 
library(msigdbr)
m_df <- msigdbr(species = "Mus musculus")
table(m_df$gs_subcat)
m_t2g <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:BIOCARTA") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "C7", subcategory = "IMMUNESIGDB") %>% 
  dplyr::select(gs_name, entrez_gene)


m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
ck <- compareCluster(geneCluster = geneid.ls, fun = enricher,TERM2GENE=m_t2g,pvalueCutoff = 0.05, 
                     #OrgDb = org.Mm.eg.db,
                     #organism = 'mmu',
                     
                     pAdjustMethod = "BH")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 
tail(ck)
#showCategory=NULL
dotplot(ck )+theme_bw()+
  theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.justification=c(1,0), 
    #legend.title = element_text("Clusters"),  
    axis.text.x  = element_text(angle=45, hjust = 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),text=element_text(family="Helvetica", size=8, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Enrichment_C19_GOBP.pdf", width = 6, height = 10)



#-----------------Cluster 1+4 : T cells + Naive T cells-----------------
clust14_TcellsNaive<-subset(x = RasVCtl_stages.combined, idents = c("T cells","Naive T cells"))
clust14_TcellsNaive

DefaultAssay(clust14_TcellsNaive) <- "integrated"
clust14_TcellsNaive <- RunPCA(clust14_TcellsNaive, verbose = T)
DimPlot(clust14_TcellsNaive, reduction = "pca")
ElbowPlot(clust14_TcellsNaive, reduction = "pca",ndims = 50)

clust14_TcellsNaive <- RunUMAP(clust14_TcellsNaive, reduction = "pca", dims = 1:50)
clust14_TcellsNaive <- FindNeighbors(clust14_TcellsNaive, reduction = "pca", dims = 1:50)
clust14_TcellsNaive <- FindClusters(clust14_TcellsNaive, resolution = 0.1)

#Assign new identities: based on both DEG and transfer guilliams
new.cluster.ids <- c("moDCs", "Classical monocytes", "BMDMs", "cDC2s","Proliferating cDC1s","Inflammatory Monocytes", "Patrolling Monocytes", "KCs",
                     "Migratory DCs")
names(new.cluster.ids) <- levels(clust14_TcellsNaive)
clust14_TcellsNaive <- RenameIdents(clust14_TcellsNaive, new.cluster.ids)


# Visualization
# Visualization
colsTcellsnaive<-sample(ArchRPalettes[[6]], size = 10, replace = F)
names(colsTcellsnaive)<-names(table(Idents(clust14_TcellsNaive)))

p1 <- DimPlot(clust14_TcellsNaive, reduction = "umap", cols = cols2[!is.na(names(cols2))], group.by = "orig.ident")#, cells = cDC1scells)+ylim(c(-15,8))+xlim(c(-5,10))#,group.by = "orig.ident")
p2 <- DimPlot(clust14_TcellsNaive, reduction = "umap",cols = colsTcellsnaive[!is.na(names(colsTcellsnaive))],label = TRUE, repel = TRUE)
p1+p2
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_controlRas12V_zoomC14_res1_annotatedclusted.pdf", width = 10, height = 5)

save(clust14_TcellsNaive, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust14_TcellsNaive.RData")
load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust14_TcellsNaive.RData")


#DEG
DefaultAssay(clust14_TcellsNaive) <- "RNA"

plan("multisession", workers = 6)
clust14_TcellsNaive.markers <- FindAllMarkers(clust14_TcellsNaive, logfc.threshold = 0.5,  min.pct =  0.25, 
                                        min.diff.pct = 0.25,assay = "RNA",slot="data")
future:::ClusterRegistry("stop")


t<-clust14_TcellsNaive.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]
c5<-t$gene[t$cluster==5]
c6<-t$gene[t$cluster==6]
c7<-t$gene[t$cluster==7]
c8<-t$gene[t$cluster==8]
c9<-t$gene[t$cluster==9]
c10<-t$gene[t$cluster==10]
c11<-t$gene[t$cluster==11]

clust14_TcellsNaive.markers_top50DEG<-cbindX( data.frame(c0), data.frame(c1),  data.frame(c2),  data.frame(c3),  data.frame(c4),  data.frame(c5),  data.frame(c6),  data.frame(c7),  data.frame(c8),  data.frame(c9))

DefaultAssay(clust14_TcellsNaive) <- "SCT"
FeaturePlot(clust14_TcellsNaive, features = c("Csf2ra","Csf2rb","Ccl6","Ccl17","Ccl24")) #C
VlnPlot(clust14_TcellsNaive, features = c("Il23a"), group.by = "orig.ident") #C
library(readxl)
Tcellsignatures_public<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/Public/Tcellsignatures.xlsx")

Tcellsignatures_public_Exhaustion<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion[!is.na(Tcellsignatures_public$Exhaustion)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_G1S<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G1/S`[!is.na(Tcellsignatures_public$`G1/S`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_G2M<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G2/M`[!is.na(Tcellsignatures_public$`G2/M`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_Exhaustion_CD8<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD8[!is.na(Tcellsignatures_public$Exhaustion_CD8)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)
Tcellsignatures_public_Exhaustion_CD4<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD4[!is.na(Tcellsignatures_public$Exhaustion_CD4)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/"), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/"), uniqueRows=T)


Tcellsignatures_public_list<-list("Memory"=Tcellsignatures_public$Memory,"Effector"=Tcellsignatures_public$Effector,
                                  "Exhaustion"=Tcellsignatures_public_Exhaustion$MGI.symbol,
                                  "G1S"=Tcellsignatures_public_G1S$MGI.symbol,  "G2M"=Tcellsignatures_public_G2M$MGI.symbol,
                                  "Exhaustion_CD8"=Tcellsignatures_public_Exhaustion_CD8$MGI.symbol,"Exhaustion_CD4"=Tcellsignatures_public_Exhaustion_CD4$MGI.symbol,
                                  "Progenitor_Exh"=Tcellsignatures_public$Progenitor_Exh,"Effector_like"=Tcellsignatures_public$Effector_like,
                                  "Terminally_Exh"=Tcellsignatures_public$Terminally_Exh, "Proliferating"=Tcellsignatures_public$Proliferating)
clust14_TcellsNaive <- AddModuleScore(
  object = clust14_TcellsNaive,
  features = Tcellsignatures_public_list,
  ctrl = 5,
  name = 'Tcellsignatures'
)
library(RColorBrewer)
p1<-FeaturePlot(object = clust14_TcellsNaive,  min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures1"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Memory")
p2<-FeaturePlot(object = clust14_TcellsNaive, min.cutoff = -5, max.cutoff = 4,features = c("Tcellsignatures2"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Effector")
p3<-FeaturePlot(object = clust14_TcellsNaive,min.cutoff = -5, max.cutoff = 4, features = c("Tcellsignatures3"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Exhaustion")
p4<-FeaturePlot(object = clust14_TcellsNaive,min.cutoff = -1, max.cutoff = 1, features = c("Tcellsignatures4"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("G1S")
p5<-FeaturePlot(object = clust14_TcellsNaive,min.cutoff = -1, max.cutoff = 1, features = c("Tcellsignatures5"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("G2M")
p6<-FeaturePlot(object = clust14_TcellsNaive,  min.cutoff = -3, max.cutoff = 3,features = c("Tcellsignatures6"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Exhaustion_CD8")
p7<-FeaturePlot(object = clust14_TcellsNaive, min.cutoff = -5, max.cutoff = 5, features = c("Tcellsignatures7"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Exhaustion_CD4")
p8<-FeaturePlot(object = clust14_TcellsNaive,min.cutoff = -5, max.cutoff = 5, features = c("Tcellsignatures8"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Progenitor_Exh")
p9<-FeaturePlot(object = clust14_TcellsNaive, min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures9"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Effector_like")
p10<-FeaturePlot(object = clust14_TcellsNaive,min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures10"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Terminally_Exh")
p11<-FeaturePlot(object = clust14_TcellsNaive, min.cutoff = -5, max.cutoff = 5,features = c("Tcellsignatures11"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("Proliferating")

pdf("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominTcellsNaiveTcells_Tcellsignatures_Enrichment.pdf", width = 13, height = 8)
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, ncol=4)
dev.off()

###########Export as excel
library(openxlsx)
SourceData <- createWorkbook()

addWorksheet(SourceData, "markerseachcluster_c14_Tcells")#Add sheets

# Write the data to the sheets
writeData(SourceData, sheet = "markerseachcluster_c14_Tcells", x = clust14_TcellsNaive.markers_top50DEG )

# Export the file
saveWorkbook(SourceData, "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominC14_TcellsNaiveTcells_DEG_MarkersClusters_RasV_res1.xlsx")






#Proportions
library(plyr)
proportions<-clust14_TcellsNaive@meta.data %>% group_by(sample_name,predicted.id, replicate) %>% summarise(orig.ident,predicted.id, replicate) %>% table()
proportions<-clust14_TcellsNaive@meta.data %>% group_by(sample_name,Idents(clust14_TcellsNaive), replicate) %>% plyr:::summarise(sample_name,Idents(clust14_TcellsNaive), replicate) %>% table()

#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:3, function(x){proportions[x,,replicate = "R1"]/t[x]})
colnames(t2)<-rownames(proportions)
#rownames(t2)<-new.cluster.ids

#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:3, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)

#R3
t5<-apply(proportions[,,replicate = "R3"], 1, sum)
t6<-sapply(1:3, function(x){proportions[x,,replicate = "R3"]/t5[x]})
colnames(t6)<-rownames(proportions)

t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)),
         data.frame("rep"="R3",reshape2::melt(t6)))
colnames(t)<-c("Rep","cluster","Model","value")


ggplot(t, aes(x=Model,y=value, fill=as.factor(cluster)) )+geom_bar(stat="identity", col="white")+theme_bw()+
  scale_fill_manual(values=cols4[!is.na(names(cols4))])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        #legend.title = element_text("Clusters"),  
        axis.text.x  = element_text(angle=0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Proportions_controlRas12D_zoomcluster19_cluster.pdf", width = 3.5, height = 6)


library(ggpubr)
library(rstatix)

t$Model<-factor(t$Model, levels = c("Ctl","RasV_Pretumor","RasV_Endpoint"))

t$type<-NA
t$type[t$cluster%in% c("Neutrophils","Monocytic cells","cDC1","pDCs","migratory DCs")]<-"Myeloid"
t$type[t$cluster%in% c("T cells","Naive T cells","B cells","NK cells","Plasma cells")]<-"Lymphoid"
t$type[t$cluster%in% c("Endothelial cells","HSC+Hepatocytes","Mixed(Hepatocytes+Fibro)","Hepatocytes")]<-"Other"

t$cluster<-as.character(t$cluster)
p <- ggbarplot(t[complete.cases(t),], x = "Model", y = "value", add = "mean_se",
               
               fill =  "cluster",palette = colsTcellsnaive[!is.na(names(colsTcellsnaive))])+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x  = element_text(angle=45, hjust=1))
p

#
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Proportions_ZoominC14_TcellsNaive_res1.pdf", width = 4, height = 6)


#Try enrichment 
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)

#Keep all the DEG
t<-clust14_TcellsNaive.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05) 
c0<-t$gene[t$avg_log2FC>0 &t$cluster==0]
c1<-t$gene[t$avg_log2FC>0 &t$cluster==1]
c2<-t$gene[t$avg_log2FC>0 &t$cluster==2]
c3<-t$gene[t$avg_log2FC>0 &t$cluster==3]
c4<-t$gene[t$avg_log2FC>0 &t$cluster==4]
c5<-t$gene[t$avg_log2FC>0 &t$cluster==5]
c6<-t$gene[t$avg_log2FC>0 &t$cluster==6]
c7<-t$gene[t$avg_log2FC>0 &t$cluster==7]
c8<-t$gene[t$avg_log2FC>0 &t$cluster==8]

clust14_TcellsNaive.markers_topUPDEGList<-list(c0,c1,  c2,  c3,c4,c5)

#Conver genename to entrez id
geneid.ls <- clust14_TcellsNaive.markers_topUPDEGList %>% map(~{
  
  
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})


names(geneid.ls)<-levels(Idents(clust14_TcellsNaive))

#Try hallmarks msigdb: 
library(msigdbr)
m_df <- msigdbr(species = "Mus musculus")
table(m_df$gs_subcat)
m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
ck <- compareCluster(geneCluster = geneid.ls, fun = enricher,TERM2GENE=m_t2g,pvalueCutoff = 0.05, 
                     #OrgDb = org.Mm.eg.db,
                     #organism = 'mmu',
                     
                     pAdjustMethod = "BH")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 
tail(ck)
#showCategory=NULL
dotplot(ck)+theme_bw()+
  theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.justification=c(1,0), 
    #legend.title = element_text("Clusters"),  
    axis.text.x  = element_text(angle=45, hjust = 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),text=element_text(family="Helvetica", size=8, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Enrichment_C19_GOBP.pdf", width = 6, height = 10)


#-----------------Cluster 3 : Monocytic cells-----------------
clust3_Monocyticcells<-subset(x = RasVCtl_stages.combined, idents = c("Monocytic cells"))
clust3_Monocyticcells

DefaultAssay(clust3_Monocyticcells) <- "integrated"
clust14_TcellsNaive <- RunPCA(clust3_Monocyticcells, verbose = T)
DimPlot(clust3_Monocyticcells, reduction = "pca")
ElbowPlot(clust3_Monocyticcells, reduction = "pca",ndims = 50)

clust3_Monocyticcells <- RunUMAP(clust3_Monocyticcells, reduction = "pca", dims = 1:50)
clust3_Monocyticcells <- FindNeighbors(clust3_Monocyticcells, reduction = "pca", dims = 1:50)
clust3_Monocyticcells <- FindClusters(clust3_Monocyticcells, resolution = 0.1)

#Assign new identities: based on both DEG and transfer guilliams
new.cluster.ids <- c("moDCs", "Classical monocytes", "BMDMs", "cDC2s","Proliferating cDC1s","Inflammatory Monocytes", "Patrolling Monocytes", "KCs",
                     "Migratory DCs")
names(new.cluster.ids) <- levels(clust3_Monocyticcells)
clust3_Monocyticcells <- RenameIdents(clust3_Monocyticcells, new.cluster.ids)

new.cluster.ids <- c("C1_Monocytic_March3+", "C2_Classical monocytes_Chil3+", "C3_Monocytic_C1qa+", "C4_cDC2_Cd209a+","C5_cDC1s_Wdfy4+", "C6_Monocytic_Cxcl3+", "C7_Monocytic_Adgre4+",
                     "C8_KC_Slc40a1+", "C9_Migratory DCs_Ccr7+")
names(new.cluster.ids) <- levels(clust3_Monocyticcells)
clust3_Monocyticcells <- RenameIdents(clust3_Monocyticcells, new.cluster.ids)

# Visualization
# Visualization
#cols_clust3_Monocyticcells<-sample(x = ArchRPalettes[[9]], size = 8, replace = F)
names(cols_clust3_Monocyticcells)<-names(table(Idents(clust3_Monocyticcells)))

p1 <- DimPlot(clust3_Monocyticcells, reduction = "umap", cols = cols2[!is.na(names(cols2))], group.by = "orig.ident")#, cells = cDC1scells)+ylim(c(-15,8))+xlim(c(-5,10))#,group.by = "orig.ident")
p2 <- DimPlot(clust3_Monocyticcells, reduction = "umap",cols = cols_clust3_Monocyticcells[!is.na(names(cols_clust3_Monocyticcells))],label = TRUE, repel = TRUE)
p1+p2
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_controlRas12V_zoomC3_Monocyticcells_res1_annotatedclusted.pdf", width = 10, height = 5)

save(clust3_Monocyticcells, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust3_Monocyticcells.RData")
load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust3_Monocyticcells.RData")

DefaultAssay(clust3_Monocyticcells) <- "SCT"
FeaturePlot(clust3_Monocyticcells, features = "Fas")
FeaturePlot(clust3_Monocyticcells, features = c("Vsig4","Clec4f","Ly6c2","Adgre1","Itgam","Folr2","Spp1"))
FeaturePlot(clust1_Tcells, features = "Fasl", reduction = "umap")
VlnPlot(clust3_Monocyticcells, features = "Pnp")


#DEG
DefaultAssay(clust3_Monocyticcells) <- "RNA"

plan("multisession", workers = 6)
clust3_Monocyticcells.markers <- FindAllMarkers(clust3_Monocyticcells, logfc.threshold = 0.5,  min.pct =  0.25, 
                                        min.diff.pct = 0.25,assay = "RNA",slot="data")
future:::ClusterRegistry("stop")


t<-clust3_Monocyticcells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]
c5<-t$gene[t$cluster==5]
c6<-t$gene[t$cluster==6]
c7<-t$gene[t$cluster==7]
c8<-t$gene[t$cluster==8]

clust3_Monocyticcells.markers_top50DEG<-cbindX( data.frame(c0), data.frame(c1),  data.frame(c2),  data.frame(c3),  data.frame(c4),  data.frame(c5),  data.frame(c6),  data.frame(c7))


###########Export as excel
library(openxlsx)
SourceData <- createWorkbook()

addWorksheet(SourceData, "markerseachcluster_c3_Mono")#Add sheets

# Write the data to the sheets
writeData(SourceData, sheet = "markerseachcluster_c3_Mono", x = clust3_Monocyticcells.markers_top50DEG )

# Export the file
saveWorkbook(SourceData, "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominC3_Monocyticcells_DEG_MarkersClusters_RasV_res1.xlsx")






#Proportions
library(plyr)
proportions<-clust3_Monocyticcells@meta.data %>% group_by(sample_name,predicted.id, replicate) %>% summarise(orig.ident,predicted.id, replicate) %>% table()
proportions<-clust3_Monocyticcells@meta.data %>% group_by(sample_name,Idents(clust3_Monocyticcells), replicate) %>% plyr:::summarise(sample_name,Idents(clust3_Monocyticcells), replicate) %>% table()

#R1
t<-apply(proportions[,,replicate = "R1"], 1, sum)
t2<-sapply(1:3, function(x){proportions[x,,replicate = "R1"]/t[x]})
colnames(t2)<-rownames(proportions)
#rownames(t2)<-new.cluster.ids

#R2
t3<-apply(proportions[,,replicate = "R2"], 1, sum)
t4<-sapply(1:3, function(x){proportions[x,,replicate = "R2"]/t3[x]})
colnames(t4)<-rownames(proportions)

#R3
t5<-apply(proportions[,,replicate = "R3"], 1, sum)
t6<-sapply(1:3, function(x){proportions[x,,replicate = "R3"]/t5[x]})
colnames(t6)<-rownames(proportions)

t<-rbind(data.frame("rep"="R1",reshape2::melt(t2)),
         data.frame("rep"="R2",reshape2::melt(t4)),
         data.frame("rep"="R3",reshape2::melt(t6)))
colnames(t)<-c("Rep","cluster","Model","value")


ggplot(t, aes(x=Model,y=value, fill=as.factor(cluster)) )+geom_bar(stat="identity", col="white")+theme_bw()+
  scale_fill_manual(values=cols4[!is.na(names(cols4))])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        #legend.title = element_text("Clusters"),  
        axis.text.x  = element_text(angle=0),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Proportions_controlRas12D_zoomcluster19_cluster.pdf", width = 3.5, height = 6)


library(ggpubr)
library(rstatix)

t$Model<-factor(t$Model, levels = c("Ctl","RasV_Pretumor","RasV_Endpoint"))

t$type<-NA
t$type[t$cluster%in% c("Neutrophils","Monocytic cells","cDC1","pDCs","migratory DCs")]<-"Myeloid"
t$type[t$cluster%in% c("T cells","Naive T cells","B cells","NK cells","Plasma cells")]<-"Lymphoid"
t$type[t$cluster%in% c("Endothelial cells","HSC+Hepatocytes","Mixed(Hepatocytes+Fibro)","Hepatocytes")]<-"Other"

t$cluster<-as.character(t$cluster)
p <- ggbarplot(t[complete.cases(t),], x = "Model", y = "value", add = "mean_se",
               
               fill =  "cluster",palette = cols_clust3_Monocyticcells[!is.na(names(cols_clust3_Monocyticcells))])+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x  = element_text(angle=45, hjust=1))
p

#
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Proportions_ZoominC3_Monocyticcells.pdf", width = 4, height = 6)


#Try enrichment 
library(clusterProfiler)
library(org.Mm.eg.db)
library(purrr)

#Keep all the DEG
t<-clust3_Monocyticcells.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj<0.05) 
c0<-t$gene[t$avg_log2FC>0 &t$cluster==0]
c1<-t$gene[t$avg_log2FC>0 &t$cluster==1]
c2<-t$gene[t$avg_log2FC>0 &t$cluster==2]
c3<-t$gene[t$avg_log2FC>0 &t$cluster==3]
c4<-t$gene[t$avg_log2FC>0 &t$cluster==4]
c5<-t$gene[t$avg_log2FC>0 &t$cluster==5]
c6<-t$gene[t$avg_log2FC>0 &t$cluster==6]
c7<-t$gene[t$avg_log2FC>0 &t$cluster==7]
c8<-t$gene[t$avg_log2FC>0 &t$cluster==8]

clust3_Monocyticcells.markers_topUPDEGList<-list(c0,c1,  c2,  c3,c4,c5,c6,c7)

#Conver genename to entrez id
geneid.ls <- clust3_Monocyticcells.markers_topUPDEGList %>% map(~{
  
  
  gene.df <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = .x,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
  
  gene <- gene.df$ENTREZID
  gene <- gene[which(!is.na(gene))]
  gene <- unique(gene)
  
  return(gene)
})


names(geneid.ls)<-levels(Idents(clust3_Monocyticcells))

#Try hallmarks msigdb: 
library(msigdbr)
m_df <- msigdbr(species = "Mus musculus")
table(m_df$gs_subcat)
m_t2g <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

m_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% 
  dplyr::select(gs_name, entrez_gene)

head(m_t2g)
ck <- compareCluster(geneCluster = geneid.ls, fun = enricher,TERM2GENE=m_t2g,pvalueCutoff = 0.05, 
                     #OrgDb = org.Mm.eg.db,
                     #organism = 'mmu',
                     
                     pAdjustMethod = "BH")
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 
tail(ck)
#showCategory=NULL
dotplot(ck)+theme_bw()+
  theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.justification=c(1,0), 
    #legend.title = element_text("Clusters"),  
    axis.text.x  = element_text(angle=45, hjust = 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),text=element_text(family="Helvetica", size=8, color="black"), axis.text =element_text(family="Helvetica", size=8, color="black"))
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Enrichment_C19_GOBP.pdf", width = 6, height = 10)




################################19. Trajectory analysis with monocle3 ############################### 
############################### 
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

DefaultAssay(clust1_Tcells)<-"SCT"
FeaturePlot(clust1_Tcells, features = "sample_name", reduction = "umap")
DimPlot(clust1_Tcells, reduction = "umap",  group.by = "sample_name", shuffle = T)
clust1_Tcells$Idents<-NULL
clust1_Tcells.cds <- as.cell_data_set(clust1_Tcells)
pData(clust1_Tcells.cds) <- NULL
#remove the control cells
#clust1_Tcells.cds <- as.cell_data_set(subset(clust1_Tcells, subset= orig.ident=="Nras12D"))
clust1_Tcells.cds <- cluster_cells(cds = clust1_Tcells.cds,resolution=1e-4, reduction_method = "UMAP")
clust1_Tcells.cds <- learn_graph(clust1_Tcells.cds, use_partition = TRUE)
plot_cells(clust1_Tcells.cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#Select as root cells
levels(Idents(clust1_Tcells))
roots<-WhichCells(subset(clust1_Tcells, subset= orig.ident=="Nras12D"),idents =  "Crip1+_Granulocytes" )
roots<-WhichCells(clust1_Tcells,idents =  "Ccl3+_Granulocytes" )

roots<-WhichCells(clust1_moMDCs,idents =  "Classical monocytes" )

roots<-WhichCells(clust19_moMDCsandDCs,idents =  "Classical monocytes" )


# order cells
clust1_Tcells.cds <- order_cells(clust1_Tcells.cds, reduction_method = "UMAP")#, root_cells = roots)
clust1_Tcells.cds <- cluster_cells(clust1_Tcells.cds)


# plot trajectories colored by pseudotime
plot_cells(
  cds = clust1_Tcells.cds,
  color_cells_by = "pseudotime",
  
  show_trajectory_graph = TRUE
)


#Try helper function for identifying roots

#It's often desirable to specify the root of the trajectory programmatically, rather than manually picking it. The function below does so by first grouping the cells according to which trajectory graph node they are nearest to. Then, it calculates what fraction of the cells at each node come from the earliest time point. Then it picks the node that is most heavily occupied by early cells and returns that as the root.

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="Crip1+_Granulocytes"){
  cell_ids <- which(colData(cds)[, "ident"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

clust1_Tcells.cds <- order_cells(clust1_Tcells.cds, root_pr_nodes=get_earliest_principal_node(clust1_Tcells.cds))
plot_cells(clust1_Tcells.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=T,
           graph_label_size=1.5)

plot_cells(clust1_Tcells.cds,
           color_cells_by = "pseudotime")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Pseudotime_Monocle3_rootCrip1bothcontrolRas12D.pdf", width = 4.5, height = 4)


get_earliest_principal_node <- function(cds, time_bin="Classical monocytes"){
  cell_ids <- which(colData(cds)[, "ident"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

clust19_moMDCsandDCs.cds <- order_cells(clust19_moMDCsandDCs.cds, root_pr_nodes=get_earliest_principal_node(clust19_moMDCsandDCs.cds))
plot_cells(clust19_moMDCsandDCs.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=T,
           graph_label_size=1.5)

plot_cells(clust1_moMDCs.cds,
           color_cells_by = "pseudotime")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Pseudotime_Monocle3_rootmoDCscontrolRas12D.pdf", width = 4.5, height = 4)

####Endothelial cluster

DefaultAssay(clust0_Endothelialcells)<-"SCT"
DimPlot(clust0_Endothelialcells, reduction = "umap",  group.by = "sample_name", shuffle = T)
clust0_Endothelialcells$Idents<-NULL
clust0_Endothelialcells.cds <- as.cell_data_set(clust0_Endothelialcells)
pData(clust0_Endothelialcells.cds) <- NULL
#remove the control cells
#clust0_Endothelialcells.cds <- as.cell_data_set(subset(clust0_Endothelialcells, subset= orig.ident=="Nras12D"))
clust0_Endothelialcells.cds <- cluster_cells(cds = clust0_Endothelialcells.cds,resolution=1e-4, reduction_method = "UMAP")
clust0_Endothelialcells.cds <- learn_graph(clust0_Endothelialcells.cds, use_partition = TRUE)
plot_cells(clust0_Endothelialcells.cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

#Select as root cells
levels(Idents(clust0_Endothelialcells))
roots<-WhichCells(subset(clust0_Endothelialcells, subset= orig.ident=="Nras12D"),idents =  "Crip1+_Granulocytes" )
roots<-WhichCells(clust0_Endothelialcells,idents =  "Ccl3+_Granulocytes" )

roots<-WhichCells(clust0_Endothelialcells,idents =  "Endothelialcells_c0" )

roots<-WhichCells(clust19_moMDCsandDCs,idents =  "Classical monocytes" )


# order cells
clust0_Endothelialcells.cds <- order_cells(clust0_Endothelialcells.cds, reduction_method = "UMAP", root_cells = roots)
clust0_Endothelialcells.cds <- cluster_cells(clust0_Endothelialcells.cds)


# plot trajectories colored by pseudotime
plot_cells(
  cds = clust0_Endothelialcells.cds,
  color_cells_by = "pseudotime",
  
  show_trajectory_graph = TRUE
)


#Try helper function for identifying roots

#It's often desirable to specify the root of the trajectory programmatically, rather than manually picking it. The function below does so by first grouping the cells according to which trajectory graph node they are nearest to. Then, it calculates what fraction of the cells at each node come from the earliest time point. Then it picks the node that is most heavily occupied by early cells and returns that as the root.

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin=""){
  cell_ids <- which(colData(cds)[, "ident"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

clust0_Endothelialcells.cds <- order_cells(clust0_Endothelialcells.cds, root_pr_nodes=get_earliest_principal_node(clust0_Endothelialcells.cds))
plot_cells(clust0_Endothelialcells.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=T,
           graph_label_size=1.5)

plot_cells(clust0_Endothelialcells.cds,
           color_cells_by = "pseudotime")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Pseudotime_Monocle3_rootCrip1bothcontrolRas12D.pdf", width = 4.5, height = 4)


get_earliest_principal_node <- function(cds, time_bin="Classical monocytes"){
  cell_ids <- which(colData(cds)[, "ident"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

clust19_moMDCsandDCs.cds <- order_cells(clust19_moMDCsandDCs.cds, root_pr_nodes=get_earliest_principal_node(clust19_moMDCsandDCs.cds))
plot_cells(clust19_moMDCsandDCs.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=T,
           graph_label_size=1.5)

plot_cells(clust1_moMDCs.cds,
           color_cells_by = "pseudotime")
ggsave("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Plots/Pseudotime_Monocle3_rootmoDCscontrolRas12D.pdf", width = 4.5, height = 4)

############################### 


##########Share with Milan lab data from endothelial cells and T cells
Endothelialcells_Tcells<-subset(x = RasVCtl_stages.combined, idents = c("Endothelial cells", "T cells"))
Endothelialcells_Tcells

DefaultAssay(Endothelialcells_Tcells) <- "integrated"
Endothelialcells_Tcells <- RunPCA(Endothelialcells_Tcells, verbose = T)
DimPlot(Endothelialcells_Tcells, reduction = "pca")
ElbowPlot(Endothelialcells_Tcells, reduction = "pca",ndims = 50)

Endothelialcells_Tcells <- RunUMAP(Endothelialcells_Tcells, reduction = "pca", dims = 1:50)
Endothelialcells_Tcells <- FindNeighbors(Endothelialcells_Tcells, reduction = "pca", dims = 1:50)
Endothelialcells_Tcells <- FindClusters(Endothelialcells_Tcells, resolution = 0.1)



# Visualization
# Visualization
cols_Endothelialcells_Tcells<-sample(x = ArchRPalettes[[3]], size = 12, replace = F)
names(cols_Endothelialcells_Tcells)<-names(table(Idents(Endothelialcells_Tcells)))
cols<-ArchRPalettes[[8]][c(3,10,11)]
names(cols)<-names(table(Endothelialcells_Tcells$sample_name))
p1 <- DimPlot(Endothelialcells_Tcells, reduction = "umap", cols = cols[!is.na(names(cols))], group.by = "sample_name", shuffle = T)#, cells = cDC1scells)+ylim(c(-15,8))+xlim(c(-5,10))#,group.by = "orig.ident")
p2 <- DimPlot(Endothelialcells_Tcells, reduction = "umap",group.by = "Idents",label = TRUE, repel = TRUE,  shuffle = T)
p1+p2
ggsave("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/UMAP_controlRas12V_Endothelialcells_Tcells.pdf", width = 10, height = 5)
Endothelialcells_Tcells$Idents<-as.character( Endothelialcells_Tcells$Idents )
save(Endothelialcells_Tcells, file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/Collab_controlRas12V_Endothelialcells_Tcells.RData")
#load(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/clust0_Endothelialcells.RData")


################## T-cell signatures enrichment ############################### 
################################
library(readxl)
Tcellsignatures_public<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/Public/Tcellsignatures.xlsx")
Tcellsignatures_Zheng_public<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/Public/Tcellsignatures_Zheng_Zhang_summarised.xlsx", sheet = 1)
Tcellsignatures_Zheng_public<-reshape2::melt(as.data.frame(Tcellsignatures_Zheng_public), measure.vars = colnames(Tcellsignatures_Zheng_public))

Tcellsignatures_Zheng_public_mouse<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values =Tcellsignatures_Zheng_public$value, mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
Tcellsignatures_Zheng_public_mouse$MGI.symbol[match(x = Tcellsignatures_Zheng_public$value, table = Tcellsignatures_Zheng_public_mouse$HGNC.symbol)]
Tcellsignatures_Zheng_public$mousegene<-Tcellsignatures_Zheng_public_mouse$MGI.symbol[match(x = Tcellsignatures_Zheng_public$value, table = Tcellsignatures_Zheng_public_mouse$HGNC.symbol)]

Tcellsignatures_public_Exhaustion<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion[!is.na(Tcellsignatures_public$Exhaustion)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105,verbose = TRUE ), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", verbose = TRUE, dataset="mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/",mirror = "asia"), uniqueRows=T)
Tcellsignatures_public_G1S<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G1/S`[!is.na(Tcellsignatures_public$`G1/S`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
Tcellsignatures_public_G2M<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$`G2/M`[!is.na(Tcellsignatures_public$`G2/M`)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
Tcellsignatures_public_Exhaustion_CD8<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD8[!is.na(Tcellsignatures_public$Exhaustion_CD8)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
Tcellsignatures_public_Exhaustion_CD4<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Tcellsignatures_public$Exhaustion_CD4[!is.na(Tcellsignatures_public$Exhaustion_CD4)], mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)


Tcellsignatures_public_list<-list("Memory"=Tcellsignatures_public$Memory,"Effector"=Tcellsignatures_public$Effector,
                                  "Exhaustion"=Tcellsignatures_public_Exhaustion$MGI.symbol,
                                  "G1S"=Tcellsignatures_public_G1S$MGI.symbol,  "G2M"=Tcellsignatures_public_G2M$MGI.symbol,
                                  "Exhaustion_CD8"=Tcellsignatures_public_Exhaustion_CD8$MGI.symbol,"Exhaustion_CD4"=Tcellsignatures_public_Exhaustion_CD4$MGI.symbol,
                                  "Progenitor_Exh"=Tcellsignatures_public$Progenitor_Exh,"Effector_like"=Tcellsignatures_public$Effector_like,
                                  "Terminally_Exh"=Tcellsignatures_public$Terminally_Exh, "Proliferating"=Tcellsignatures_public$Proliferating)

Tcellsignatures_Zheng_list<-list("C1_CD8−LEF1"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C1_CD8−LEF1"],
                                 "C2_CD8-CX3CR1"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C2_CD8-CX3CR1"],
                                 "C3_CD8-SLC4A10"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C3_CD8-SLC4A10"],
                                 "C4_CD8-LAYN"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C4_CD8-LAYN"],
                                 "C5_CD8-GZMK"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C5_CD8-GZMK"],
                                 "C6_CD4-CCR7"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C6_CD4-CCR7"],
                                 "C7_CD4-FOXP3"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C7_CD4-FOXP3"],
                                 "C8_CD4-CTLA4"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C8_CD4-CTLA4"],
                                 "C9_CD4-GZMK"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C9_CD4-GZMK"],
                                 "C10_CD4-CXCL13"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C10_CD4-CXCL13"],
                                 "C11_CD4-GNLY"=Tcellsignatures_Zheng_public$mousegene[!is.na(Tcellsignatures_Zheng_public$mousegene) & Tcellsignatures_Zheng_public$variable=="C11_CD4-GNLY"])
clust0_Endothelialcells <- AddModuleScore(
  object = clust1_Tcells,
  features = Tcellsignatures_Zheng_list,
  ctrl = 5,
  name = 'Tcellsignatures_Zheng'
)

library(RColorBrewer)
p1<-FeaturePlot(object = clust1_Tcells,  min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures_Zheng1"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C1_CD8−LEF1")#+ggtitle("Memory")
p2<-FeaturePlot(object = clust1_Tcells, min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures_Zheng2"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C2_CD8-CX3CR1")#+ggtitle("Effector")
p3<-FeaturePlot(object = clust1_Tcells,min.cutoff = -2, max.cutoff = 2, features = c("Tcellsignatures_Zheng3"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C3_CD8-SLC4A10")#+ggtitle("Exhaustion")
p4<-FeaturePlot(object = clust1_Tcells,min.cutoff = -2, max.cutoff = 1.5, features = c("Tcellsignatures_Zheng4"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C4_CD8-LAYN")#+ggtitle("G1S")
p5<-FeaturePlot(object = clust1_Tcells,min.cutoff = -2, max.cutoff = 2, features = c("Tcellsignatures_Zheng5"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C5_CD8-GZMK")#+ggtitle("G2M")
p6<-FeaturePlot(object = clust1_Tcells,  min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures_Zheng6"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C6_CD4-CCR7")#+ggtitle("Exhaustion_CD8")
p7<-FeaturePlot(object = clust1_Tcells, min.cutoff = -2, max.cutoff = 2, features = c("Tcellsignatures_Zheng7"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C7_CD4-FOXP3")#+ggtitle("Exhaustion_CD4")
p8<-FeaturePlot(object = clust1_Tcells,min.cutoff = -2, max.cutoff = 2, features = c("Tcellsignatures_Zheng8"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C8_CD4-CTLA4")#+ggtitle("Progenitor_Exh")
p9<-FeaturePlot(object = clust1_Tcells, min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures_Zheng9"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C9_CD4-GZMK")#+ggtitle("Effector_like")
p10<-FeaturePlot(object = clust1_Tcells,min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures_Zheng10"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C10_CD4-CXCL13")
p11<-FeaturePlot(object = clust1_Tcells, min.cutoff = -2, max.cutoff = 2,features = c("Tcellsignatures_Zheng11"))+     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+ggtitle("C11_CD4-GNLY")

pdf("/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/ZoominTcells_TcellsignaturesZheng_Enrichment.pdf", width = 13, height = 8)
gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, ncol=4)
dev.off()
################################


################################LAMs signature enrichment ###############################
################################
library(readxl)
library(hypeR)

AT_LAMs<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Public/Guilliams2022/LAMsSignature.xlsx",sheet = 1,col_names = F)
AT_LAMs$...1

HepvsAT_LAMs<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Public/Guilliams2022/LAMsSignature.xlsx",sheet = 2,col_names = T)

HepvsAT_LAMs %>% filter(abs(pct.1-pct.2)>0.5 &avg_logFC > 0 &  p_val_adj<=0.05) %>% summarise(gene) ->Hep_LAMs

LAMs_liver<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Public/Guilliams2022/LAMsSignature_Guilliams.xlsx",sheet = 5,col_names = T)

genesets<-list("AT_LAMs"=AT_LAMs[[1]], "Hep_LAMs"=Hep_LAMs[,1],"Liver_LAMs"=LAMs_liver[[1]] )
#DEG per cluster significant and logFC>0
t<-clust3_Monocyticcells.markers %>%
  group_by(cluster) %>% filter(p_val_adj<=0.5 & avg_log2FC>0) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]
c5<-t$gene[t$cluster==5]
c6<-t$gene[t$cluster==6]
c7<-t$gene[t$cluster==7]
clust3_Monocyticcells.markers_top50DEG_Upregulated<-cbindX(as.data.frame(c0),as.data.frame(c1),as.data.frame(c2),as.data.frame(c3),as.data.frame(c4),as.data.frame(c5),as.data.frame(c6),as.data.frame(c7))


clust3_Monocyticcells.markers_top50DEG_Upregulated<-list((c0),(c1),(c2),(c3),(c4),(c5),(c6),(c7))
clust3_Monocyticcells.markers_top50DEG_Upregulated
names(clust3_Monocyticcells.markers_top50DEG_Upregulated)<-c("c0","c1","c2","c3","c4","c5","c6","c7")
hyp_obj <- hypeR(signature = clust3_Monocyticcells.markers_top50DEG_Upregulated, genesets = genesets, test = "hypergeometric")
hyp_dots(hyp_obj,merge=TRUE,sizes = T,fdr = 0.05,title="LAMs")+
  #scale_color_distiller( palette = "Greens", direction = -1)+
  theme(text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
hyp_Monocytic.combined<-hyp_obj$as.list()
hyp_Monocytic.combined2<-do.call(rbind.data.frame, lapply(hyp_Monocytic.combined, function(x){x[c(1:3,6)]}))
hyp_Monocytic.combined2$cluster<-gsub(pattern = ".AT_LAMs|.Liver_LAMs|.Hep_LAMs",replacement = "", rownames(hyp_Monocytic.combined2))
ggplot(data = hyp_Monocytic.combined2, aes(x=reorder(cluster, -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  facet_grid(~label,  scales = "free")+
  #  scale_y_discrete(name="",breaks=c("Myc_p53minusControl","Myc_PTENminusControl","Ras12D_PTENminusControl","Ras12V_PTENminusControl"),
  #                   labels=c(expression("Myc" ^"OE"~"Trp53" ^"KO"), expression("Myc" ^"OE"~"Pten" ^"KO"),expression("Nras12D" ^"OE"~"\nPten" ^"KO"),expression("Nras12V" ^"OE"~"\nPten" ^"KO")))+
  ggtitle("LAMs")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Reds", direction = 1)+
  #scale_fill_viridis(option = "B",begin = .1, end = .9,direction = -1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/LAMs_Monocytic.pdf", width = 6, height = 3)



################################



################################PLVAP ECs signature enrichment ###############################
################################
library(readxl)
library(hypeR)
oncofetal_ECs<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Public/Sharma2020/DEGpercluster.xlsx",sheet = 2,col_names = T, skip = 2)
oncofetal_ECs$`Cluster 4`

oncofetal_ECs_public_mouse<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values =oncofetal_ECs$`Cluster 4`, mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
oncofetal_ECs_public_mouse$MGI.symbol[match(x = oncofetal_ECs$`Cluster 4`, table = oncofetal_ECs_public_mouse$HGNC.symbol)]
oncofetal_ECs_public_mouse$mousegene<-oncofetal_ECs_public_mouse$MGI.symbol[match(x = oncofetal_ECs$`Cluster 4`, table = oncofetal_ECs_public_mouse$HGNC.symbol)]


genesets<-list("PLVAP+ECs"=oncofetal_ECs_public_mouse$MGI.symbol)
#DEG per cluster significant and logFC>0
t<-clust0_Endothelialcells.markers %>%
  group_by(cluster) %>% filter(p_val_adj<=0.5 & avg_log2FC>0) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]

clust0_Endothelialcells.markers_top50DEG_Upregulated<-cbindX(as.data.frame(c0),as.data.frame(c1),as.data.frame(c2),as.data.frame(c3),as.data.frame(c4))


clust0_Endothelialcells.markers_top50DEG_Upregulated<-list((c0),(c1),(c2),(c3),(c4))
clust0_Endothelialcells.markers_top50DEG_Upregulated
names(clust0_Endothelialcells.markers_top50DEG_Upregulated)<-c("c0","c1","c2","c3","c4")
hyp_obj <- hypeR(signature = clust0_Endothelialcells.markers_top50DEG_Upregulated, genesets = genesets, test = "hypergeometric")
hyp_dots(hyp_obj,merge=TRUE,sizes = T,fdr = 0.05)+
  #scale_color_distiller( palette = "Greens", direction = -1)+
  theme(text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
hyp_Endothelialcells.combined<-hyp_obj$as.list()
hyp_Endothelialcells.combined2<-do.call(rbind.data.frame, lapply(hyp_Endothelialcells.combined, function(x){x[c(1:3,6)]}))
hyp_Endothelialcells.combined2$cluster<-rownames(hyp_Endothelialcells.combined2)
ggplot(data = hyp_Endothelialcells.combined2, aes(x=reorder(cluster, -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  facet_grid(~label,  scales = "free")+
  #  scale_y_discrete(name="",breaks=c("Myc_p53minusControl","Myc_PTENminusControl","Ras12D_PTENminusControl","Ras12V_PTENminusControl"),
  #                   labels=c(expression("Myc" ^"OE"~"Trp53" ^"KO"), expression("Myc" ^"OE"~"Pten" ^"KO"),expression("Nras12D" ^"OE"~"\nPten" ^"KO"),expression("Nras12V" ^"OE"~"\nPten" ^"KO")))+
  ggtitle("PLVAP+ECs")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Reds", direction = 1)+
  #scale_fill_viridis(option = "B",begin = .1, end = .9,direction = -1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/PLVAP_Endothelialcells.pdf", width = 3, height = 3)



################################



################################ ECs signature enrichment ###############################
################################
library(readxl)
library(hypeR)
TECs<-read_excel("/Users/m.ando/surfdrive/Documents/HCC/scRNAseq/Public/Goveia2020/degpercluster.xlsx",sheet = 1,col_names = T)

TECs_public_mouse<-getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values =TECs$Feature, mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version = 105), attributesL = c("mgi_symbol"), martL = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version = 105), uniqueRows=T)
TECs_public_mouse$MGI.symbol[match(x = TECs$Feature, table = TECs_public_mouse$HGNC.symbol)]
TECs$mousegene<-TECs_public_mouse$MGI.symbol[match(x = TECs$Feature, table = TECs_public_mouse$HGNC.symbol)]
names(table(TECs$Type))

genesets<-list("Activated capillary"=TECs$mousegene[TECs$Type=="activated capillary" & !is.na(TECs$mousegene)],"Activated pcv"=TECs$mousegene[TECs$Type=="activated pcv"& !is.na(TECs$mousegene)],
               "Artery"=TECs$mousegene[TECs$Type=="artery"& !is.na(TECs$mousegene)],"Capillary type I"=TECs$mousegene[TECs$Type=="capillary type I"& !is.na(TECs$mousegene)],
               "Capillary type II"=TECs$mousegene[TECs$Type=="capillary type II"& !is.na(TECs$mousegene)],"Immature"=TECs$mousegene[TECs$Type=="immature"& !is.na(TECs$mousegene)],
               "Postcapillary vein"=TECs$mousegene[TECs$Type=="postcapillary vein"& !is.na(TECs$mousegene)],"Scavenging capillary"=TECs$mousegene[TECs$Type=="scavenging capillary"& !is.na(TECs$mousegene)],
               "Tip cell"=TECs$mousegene[TECs$Type=="tip cell"& !is.na(TECs$mousegene)])
#DEG per cluster significant and logFC>0
t<-clust0_Endothelialcells.markers %>%
  group_by(cluster) %>% filter(p_val_adj<=0.5 & avg_log2FC>0) %>%
  slice_max(n = 50, order_by = avg_log2FC)

c0<-t$gene[t$cluster==0]
c1<-t$gene[t$cluster==1]
c2<-t$gene[t$cluster==2]
c3<-t$gene[t$cluster==3]
c4<-t$gene[t$cluster==4]

clust0_Endothelialcells.markers_top50DEG_Upregulated<-cbindX(as.data.frame(c0),as.data.frame(c1),as.data.frame(c2),as.data.frame(c3),as.data.frame(c4))


clust0_Endothelialcells.markers_top50DEG_Upregulated<-list((c0),(c1),(c2),(c3),(c4))
clust0_Endothelialcells.markers_top50DEG_Upregulated
names(clust0_Endothelialcells.markers_top50DEG_Upregulated)<-c("c0","c1","c2","c3","c4")
hyp_obj <- hypeR(signature = clust0_Endothelialcells.markers_top50DEG_Upregulated, genesets = genesets, test = "hypergeometric")
hyp_dots(hyp_obj,merge=TRUE,sizes = T,fdr = 0.05)+
  #scale_color_distiller( palette = "Greens", direction = -1)+
  theme(text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
hyp_Endothelialcells.combined<-hyp_obj$as.list()
hyp_Endothelialcells.combined2<-do.call(rbind.data.frame, lapply(hyp_Endothelialcells.combined, function(x){x[c(1:3,6)]}))
hyp_Endothelialcells.combined2$cluster<-rownames(hyp_Endothelialcells.combined2)
ggplot(data = hyp_Endothelialcells.combined2, aes(x=reorder(cluster, -1*log10(fdr)), y=-1*log10(fdr),fill=overlap))+geom_bar(stat="identity", col="black")+coord_flip()+theme_bw()+
  facet_grid(~label,  scales = "free")+
  #  scale_y_discrete(name="",breaks=c("Myc_p53minusControl","Myc_PTENminusControl","Ras12D_PTENminusControl","Ras12V_PTENminusControl"),
  #                   labels=c(expression("Myc" ^"OE"~"Trp53" ^"KO"), expression("Myc" ^"OE"~"Pten" ^"KO"),expression("Nras12D" ^"OE"~"\nPten" ^"KO"),expression("Nras12V" ^"OE"~"\nPten" ^"KO")))+
  ggtitle("PLVAP+ECs")+
  geom_hline(yintercept = 1.3)+
  scale_fill_distiller( palette = "Reds", direction = 1)+
  #scale_fill_viridis(option = "B",begin = .1, end = .9,direction = -1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification=c(1,0), 
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.background = element_blank(),text=element_text(family="Helvetica", size=12, color="black"), axis.text =element_text(family="Helvetica", size=12, color="black"))
ggsave(file = "/Users/m.ando/surfdrive/Documents/HCC/NASH/scRNAseq/CANIL/PLVAP_Endothelialcells.pdf", width = 3, height = 3)



################################
