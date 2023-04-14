##    title: "SCRGOT 2023 Coder Upgrade Session 4: Comining single cell data sets"
##    author: "Xin Wang"
##    date: '2023-5-3'
##    email: xin.wang@nationwidechildrens.org
##    output: html_document

## Method Part1: Integrating using the Harmony 
## Method reference: https://satijalab.org/seurat/articles/integration_introduction.html

## Dataset Description:
## This tutorial has four samples, 
## including two replicates from human normal kidney (GSM5627690, GSM5627691), 
## two replicates from  autosomal dominant polycystic kidney disease (GSM5627695 and GSM5627696).
## The datasets can be downloaded from GEO:GSE185948 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185948
## To reduce the time, We randomly select 1000 of cells in each sample.
## Dataset reference: Muto Y, Dixon EE, Yoshimura Y, Wu H et al. Defining cellular complexity in human autosomal dominant polycystic kidney disease by multimodal single cell analysis. Nat Commun 2022 Oct 30;13(1):6497. PMID: 36310237

## please download the datasets from the github 
## github datasets: 


############################################
# removing batche effects with Harmony
############################################

## Step 1: loading the packages, 
## Please install all the packages before loading them.

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

## set up the workspace environment:
setwd("/Users/XXW004/Documents/Projects/LectureForSingleCellGroup/IntegrationSection/TestDatasets/")

## check the files in the current directory
Datafiles<-list.files(path = "./", recursive = F, full.names = F)

## 
## Using loop to read the objects, add group column for each meta object, and save the object into a list accordingly to their names.
## 
## Create an empty list to store the objects
myobjectlist<- list()

for (x in Datafiles){
  ## find the file with the substring of ".small.rds", then replace with empty string
  #x=c("GSM5627695_PKD1.small.rds")
  name<-gsub (".small.rds","", x)
  ## loading the RDS file
  rds<- readRDS(file = x)
  dim(rds)
  head(rds@meta.data)
  # length(rds$orig.ident)
  ## using the assign the name to the objects
  assign (name, rds)
  ## create a group that show the name of samples
  rds$group <- rep(name, length(rds$orig.ident))
  myobjectlist[[name]] <-rds
}

# check the lists 
myobjectlist

# Check how many objects in the list
length(myobjectlist)
# Check the first object meta data
myobjectlist[[1]]@meta.data

## Merge multiple objects from the list, add cell id with different prefix, and create a project for the merged object.
scrna<-merge(x=myobjectlist[[1]], y=c(myobjectlist[[2]], myobjectlist[[3]],myobjectlist[[4]]), add.cell.ids = c("A","B","C","D"), project="Integration")

# Check the structure of meta data
str(scrna@meta.data)

# View the meta data
View(scrna@meta.data)

## QC & filtering
# calculate mitochondrial percentatge

scrna$mitoPercent <-PercentageFeatureSet(scrna, pattern = '^MT-')

# let's check the quality of datasets by eveluating the mitochondrial percentage, number of Counts, number of features.
head(scrna@meta.data)
VlnPlot(scrna, features = c("mitoPercent", "nCount_RNA", "nFeature_RNA"))
VlnPlot(scrna, features = c("mitoPercent", "nCount_RNA", "nFeature_RNA") , split.by = 'group')

# filtering

scrna <- subset (scrna, subset =mitoPercent <10 & nFeature_RNA >500 & nCount_RNA >200 )

# # split the object
# InstallData("ifnb")
# 
# # load dataset
# LoadData("ifnb")
# 
# ifnb@meta.data
# # split the objects


# perform the standard workflow to figure out if there are any batch effects
scrna<- NormalizeData(object = scrna)

scrna<- FindVariableFeatures(object = scrna, nfeatures = 3000)

scrna<- ScaleData(object = scrna)

scrna<- RunPCA(object = scrna, npcs =15)

scrna<- RunUMAP(scrna, dims = 1:15)
BeforeHarmony<- DimPlot(scrna, reduction = "umap",split.by = 'group')
BeforeHarmony
## run harmony

scrna@meta.data
seurat.harmony.integrated <- RunHarmony(scrna, group.by.vars = 'group', dims.use = 1:15, plot_convergence= FALSE, project.dim = F)

seurat.harmony.integrated@reductions

seurat.harmony.integrated.embed <- Embeddings(seurat.harmony.integrated, "harmony")

seurat.harmony.integrated.embed[1:10,1:10]

# Run Umap and clusering using Harmony reduction
?RunUMAP
seurat.harmony.integrated <- RunUMAP(seurat.harmony.integrated,reduction='harmony', dim=1:15)

seurat.harmony.integrated<-FindNeighbors(seurat.harmony.integrated, reduction = 'harmony', dims = 1:15)

seurat.harmony.integrated<- FindClusters(seurat.harmony.integrated, resolution=1)

# visualization

DimHarmony<-DimPlot(seurat.harmony.integrated, reduction = 'umap', split.by = 'group')
DimHarmony
## let's also compare with the CCA method
# find integration anchor (CCA)
## for four splited objects to perform the CCA method
# split the objects
SplitedObjects<- SplitObject(scrna, split.by = 'group')

# check the split objects
SplitedObjects

length(SplitedObjects)

## Normalized dataset and Find variable features
for (i in 1: length(SplitedObjects)){
  SplitedObjects[[i]] <-NormalizeData(object = SplitedObjects[[i]])
  SplitedObjects[[i]] <- FindVariableFeatures(object = SplitedObjects[[i]],selection.method = "vst")
}


# select  intergration features

features<- SelectIntegrationFeatures(object.list = SplitedObjects)
head(features)

# find integration anchor (CCA)

anchors<- FindIntegrationAnchors(object.list = SplitedObjects, anchor.features = features)

#anchors<- FindIntegrationAnchors(object.list = SplitedObjects, anchor.features = features)
# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)
?IntegrateData
# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
# we can see from the PCA that a good overlay of several condtions by PCA
seurat.integrated<- FindNeighbors(seurat.integrated, dims = 1:15, reduction = "pca")
seurat.integrated <- FindClusters(seurat.integrated)
# Now, we can also visualize with UMAP.

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15, reduction = "pca")

# Plot UMAMP
DimPlot(seurat.integrated)

#Idents(seurat.integrated)<-ordered(factor(Idents(seurat.integrated)), levels=c(0:16))

DimCCA <- DimPlot(seurat.integrated, reduction = 'umap', split.by  = 'group')

## compare two methods
grid.arrange(DimHarmony, DimCCA,  ncol = 1, nrow=2)


## potential way to check the integration methods:
## using the biomarkers

FeaturePlot(seurat.harmony.integrated, features = c("PECAM1"), split.by = 'group', min.cutoff = 0.1)

FeaturePlot(seurat.integrated, features = c("PECAM1"), split.by = 'group', min.cutoff = 0.1)

?FeaturePlot
