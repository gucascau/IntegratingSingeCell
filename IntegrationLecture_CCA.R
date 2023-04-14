##    title: "SCRGOT 2023 Coder Upgrade Session 4: Comining single cell data sets"
##  author: "Xin Wang"
##  date: '2023-5-3'
##     email: xin.wang@nationwidechildrens.org
##   output: html_document

## Method Part1: Integrating using the canonical correlation analysis (CCA) in Seruat
## Method reference: https://satijalab.org/seurat/articles/integration_introduction.html

## Dataset Description:
## This tutorial has four samples, 
## including two replicates from human normal kidney (GSM5627690, GSM5627691), 
## two replicates from  autosomal dominant polycystic kidney disease (GSM5627695 and GSM5627696).
## The dataset can be downloaded from GEO:GSE185948 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185948
## To reduce the time, We randomly select 1000 of cells in each sample.
## Dataset reference: Muto Y, Dixon EE, Yoshimura Y, Wu H et al. Defining cellular complexity in human autosomal dominant polycystic kidney disease by multimodal single cell analysis. Nat Commun 2022 Oct 30;13(1):6497. PMID: 36310237

## please download the datasets from the github 
## github datasets: 


## Step 1: loading the packages, 
## Please install all the packages before loading them.
library(Seurat)
library(gridExtra)
library(ggplot2)
library(tidyverse)

## Step 2: set up the workspace environment:
setwd("/Users/XXW004/Documents/Projects/LectureForSingleCellGroup/IntegrationSection/TestDatasets/")

## check the files in the current directory
Datafiles<-list.files(path = "./", recursive = F, full.names = F)
Datafiles

## Step 3: reading each sample

Control1<- readRDS(file = "GSM5627690_cont1.small.rds")
Control2<- readRDS(file = "GSM5627691_cont2.small.rds")
Disease1<- readRDS(file = "GSM5627695_PKD1.small.rds")
Disease2<- readRDS(file = "GSM5627696_PKD2.small.rds")

## briefly check wether we successfully load the single cell files
dim(Control1@meta.data)
str(Control1@meta.data)
head(Control1@meta.data)


## merge multiple objects from the list, we also set up the cell id 
scrna<-merge(x=Control1, y=c(Control2, Disease1,Disease2), add.cell.ids = c("A","B","C","D"), project="Integration")
scrna@meta.data


## we then create a group column fro the meta data based on their condtion

## mutate to add a group into the meta, 
## str_split split the cell ID and select the first element, simplify =TRUE will return a character matrix

scrna@meta.data <- scrna@meta.data %>% mutate(group = stringr::str_split(row.names(scrna@meta.data), "_", simplify = TRUE)[,1])

## please check the string split function:
test<-str_split(row.names(scrna@meta.data), "_", simplify = TRUE)[,1]
tail(scrna@meta.data)
# check the structure of meta data
str(scrna@meta.data)
levels(factor(scrna@meta.data$group))
# view the meta data
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

# perform the standard workflow to figure out if there are any batch effects
scrna<- NormalizeData(object = scrna)

scrna<- FindVariableFeatures(object = scrna)

scrna<- ScaleData(object = scrna)

scrna<- RunPCA(object = scrna)

ElbowPlot(scrna)

scrna<- FindNeighbors(object = scrna, dim=1:15)
scrna<- FindClusters(object=scrna)
scrna<- RunUMAP(object = scrna, dims = 1:15)


# plot 
p1 <- DimPlot(scrna, reduction = 'umap', group.by = 'group')
p1

p2<- DimPlot(scrna, reduction = 'umap', split.by = 'group')

grid.arrange(p1, p2, ncol=2, nrow=1)

# Obviously, there are batch effects for the sample, so we perform integration to correct for batch effects

# split the objects
SplitedObjects<- SplitObject(scrna, split.by = 'group')

# check the split objects
SplitedObjects

length(SplitedObjects)

## Normalized dataset and Find variable features
for (i in 1: length(SplitedObjects)){
  SplitedObjects[[i]] <-NormalizeData(object = SplitedObjects[[i]])
  SplitedObjects[[i]] <- FindVariableFeatures(object = SplitedObjects[[i]])
}

# select  intergration features

features<- SelectIntegrationFeatures(object.list = SplitedObjects)
head(features)

# find integration anchor (CCA)

anchors<- FindIntegrationAnchors(object.list = SplitedObjects, anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
# Plot PCA
PCAPlot(seurat.integrated, split.by= "group")

# we can see from the PCA that a good overlay of several condtions by PCA

# Now, we can also visualize with UMAP.

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15, reduction = "pca")

# Plot UMAMP
DimPlot(seurat.integrated)

#Idents(seurat.integrated)<-ordered(factor(Idents(seurat.integrated)), levels=c(0:16))

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'group')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'group',
              cols = c('red','green','blue','yellow'))
#p5 <- DimPlot(seurat.integrated, reduction = 'umap', split.by = 'group',
 #             cols = c('red','green','blue','yellow'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow=2)

## side by side splitted by groups

DimPlot(seurat.integrated, reduction = 'umap', split.by = 'group')


############################################
# removing batche effects with Harmony
############################################

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)

# after the merge

Datafiles<-list.files(path = "./", recursive = F, full.names = F)

Datafiles
# ### read the h5 files
# 
# for (x in Datafiles){
#   #x= "GSM5627690_cont1_filtered_feature_bc_matrix.h5"
#   # get the name file
#   name<- gsub("_filtered_feature_bc_matrix.h5","",x)
#   name
#   singlecell<-Read10X_h5(x)
#   
#   head(singlecell)
#   # create an seurate object
#   seurat_ojects<- CreateSeuratObject(counts = singlecell)
#   seurat_ojects
#   # for reproducibiltiy, set a random seed
#   set.seed(1111)
#   
#   
#   # random subset the object
#   sample_random<- subset(seurat_ojects, cells = sample(Cells(seurat_ojects),1000))
# 
#   
#   
#   
#   saveRDS(sample_random, file = paste0(name,".small.rds"))
#   
# }

Datafiles


## using the for loop to read the objects
## We can save these objects in the list
## create an empty list
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
  ## add a group that show the name of samples
  rds$group <- rep(name, length(rds$orig.ident))
  myobjectlist[[name]] <-rds
  
}

# check the lists 
myobjectlist

# how many objects in the list
length(myobjectlist)
# check the first object
myobjectlist[[1]]@meta.data


## merge multiple objects from the list
scrna<-merge(x=myobjectlist[[1]], y=c(myobjectlist[[2]], myobjectlist[[3]],myobjectlist[[4]]), add.cell.ids = c("A","B","C","D"), project="Integration")


# check the structure of meta data
str(scrna@meta.data)

# view the meta data
View(scrna@meta.data)

# 


## QC & filtering

# calculate mitochondrial percentatge

scrna$mitoPercent <-PercentageFeatureSet(scrna, pattern = '^MT-')

# let's check the quality of datasets by eveluating the mitochondrial percentage, number of Counts, number of features.
head(scrna@meta.data)
VlnPlot(scrna, features = c("mitoPercent", "nCount_RNA", "nFeature_RNA"))
VlnPlot(scrna, features = c("mitoPercent", "nCount_RNA", "nFeature_RNA") , split.by = 'group')

# filtering

scrna <- subset (scrna, subset =mitoPercent <10 & nFeature_RNA >500 & nCount_RNA >200 )

# perform the standard workflow to figure out if there are any batch effects
scrna<- NormalizeData(object = scrna)

scrna<- FindVariableFeatures(object = scrna)

scrna<- ScaleData(object = scrna)

scrna<- RunPCA(object = scrna)
scrna@meta.data
# run harmony
seurat.harmony.integrated <- RunHarmony(scrna, group.by.vars = 'group', plot_convergence= FALSE)

seurat.harmony.integrated@reductions

seurat.harmony.integrated.embed <- Embeddings(seurat.harmony.integrated, "harmony")

seurat.harmony.integrated.embed [1:40,1:40]

# Run Umap and clusering using Harmony reduction
?RunUMAP
seurat.harmony.integrated <- RunUMAP(seurat.harmony.integrated,reduction='harmony', dim=1:15)

seurat.harmony.integrated<-FindNeighbors(seurat.harmony.integrated, reduction = 'harmony', dims = 1:15)

seurat.harmony.integrated<- FindClusters(seurat.harmony.integrated, resolution=1)

# visualization

DimPlot(seurat.harmony.integrated, reduction = 'umap', split.by = 'group')

# 
