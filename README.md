# Xenium Spatial Analysis R Toolkit  
<img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/e3f99642-03a7-4b90-ae9e-a85eab71bb97" width="20%"></img>    
This repository contains R tools and workflows to spatially analyze cell-types using single-cell spatial transcriptomic *Xenium* data.   

## Table of Contents
1. [Prep Environment](#prep-environment-installing-seurat)
2. [Tutorial](#tutorial)
3. [Toolkit Contents](#toolkit-contents)
4. [Function List](#function-list)


# Prep Environment: installing Seurat
These tools are compatible with Seurat V5. For more information, see https://satijalab.org/seurat/articles/install_v5.   
To install and load package:
```R
install.packages('Seurat')
library(Seurat)
```

# Tutorial 
 
## 1. Read in Data ##
#### Read in raw Xenium data #### 
Read in raw Xenium data using Seurat's `LoadXenium()` function. Make sure to set unique FOV names for each slice for easier downstream processing. The code here also remoces cells with zero counts, and saves data as .rds file in a specified directory. 

For each replicate assign a unique FOV for easier downstream processing:
```R
#Margarita Kapustina, UBC 2023
#read in xenium data, and remove cells with 0 counts, save rds files.

#load packages
library(Seurat)
library(dplyr)

#specify folder location (home) and where to save .rds object (save.location)
home = "XeniumOutputFolder.path"
save.location = "SaveLocation.path"

#for each replicate assign a unique FOV for easier downstream processing:
setwd(home)
filename = 'output-XETG00098__0018569__rep1_atn_ant__20231129__015936'
obj <- LoadXenium(filename, fov = "fov")
# remove cells with 0 counts
obj <- subset(obj, subset = nCount_Xenium > 0)
setwd(save.location)
saveRDS(obj, 'output-XETG00098__0018569__rep1_atn_ant__20231129__015936-obj.rds')
rm(obj)
```
## 2. Merge data ##
#### Merge data across slices and animals #### 
Read in data from .rds files, add sample metadata to differentiate slices, and merge slices. 
```R
#load libraries
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

#read in data from saved .rds files and merge all slices  
data.dir = "X:/Cembrowski Lab/Margo/ATN-23/xenium-rds-objects-uniqueFOV/output-XETG00098__0018569__rep1_atn_ant__20231129__015936-obj.rds"
S1_ant = readRDS(data.dir)
data.dir = "X:/Cembrowski Lab/Margo/ATN-23/xenium-rds-objects-uniqueFOV/output-XETG00098__0018569__rep1_atn_int__20231129__015936-obj.rds"
S1_int = readRDS(data.dir)
data.dir = "X:/Cembrowski Lab/Margo/ATN-23/xenium-rds-objects-uniqueFOV/output-XETG00098__0018569__rep1_atn_post__20231129__015936-obj.rds"
S1_post = readRDS(data.dir)

#add sample metadata
S1_ant@meta.data$sample = 'S1_ant'
S1_int@meta.data$sample = 'S1_int'
S1_post@meta.data$sample = 'S1_post'

#merge slices
merge = merge(S1_ant, S1_int)
merge = merge(merge, S1_post)
merge = merge(merge, S2_ant)
```

## 3. Analyze data ##
#### 3a. Analyze merged data #### 
Analyze data across multiple slices, via [UMAP dimensionality reduction](https://www.nature.com/articles/nbt.4314) and applying [cluster-based algorithms](https://www.tandfonline.com/doi/full/10.1080/15476286.2020.1728961). Continue until you have desired UMAP embeddings. 
For more info see [Seurat's vigentte](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2.html).

```R
#normalize merged data
merge <- SCTransform(merge, assay = "Xenium")
#remove oligos and endothelial cells - your genes of choice here
xenexc = merge %>% subset(Slc17a6 > 0 & Opalin ==0 & Sox10 ==0 & Sox17 ==0)
#vizualize subset of cells in specified FOV
ImageDimPlot(xenexc, fov = 'X1fov')

#run analysis
xenexc <- RunPCA(xenexc, npcs = 30, features = rownames(xenexc))
#select your number of dims based on elbow plot
     #fyi elbow plot: a ranking of principle components based on the percentage of variance explained by each one (see: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
ElbowPlot(xenexc, ndims = 30, reduction = "pca") 
xenexc <- RunUMAP(xenexc, dims = 1:30)
xenexc <- FindNeighbors(xenexc, reduction = "pca", dims = 1:30)
xenexc <- FindClusters(xenexc, resolution = 0.03)

#visualize UMAP
DimPlot(xenexc, raster = FALSE, cols = 'glasbey') + DarkTheme() +coord_fixed()
```
To highlight a cluster in your Xenium object, use `HighlightCluster()`.
```R
highlightCluster(xenexc, cluster_id = 3)
**embed!!***
```
#### 3a. Subset clusters of choice for downstream analysis #### 
Subset clusters from a Xenium object with `subset()`.
```R
xen_atn = xenexc %>% subset(idents = c("2", "6"))
```
To highlight your subset object within the original object, use `HighlightCells()`.
```R
highlightCluster(highlight_obj = xenatn, within_obj = xenexc)
```
> [!TIP]
> #### Repeat these steps until you are satisfied with your final subset of cells. ####
(ie example ATN cells!)

## 4. Begin spatial analysis! ##
#### Create marker gene boxplots for each cluster #### 
Create boxplots for gene of interest across all clusters, using `XeniumBoxPlot()` or `XeniumBoxPlotRaw()`
```R
#create vector of genes
my_genes= c('C1ql2','Slc17a7', 'Gng13')

#generat boxplots using sequencing depth corrected counts (SCT$counts)
BoxPlots = XeniumBoxPlot(object = atn, genes = my_genes)

#or using raw counts (Xenium$counts)
BoxPlotsRaw = XeniumBoxPlot(object = atn, genes = my_genes)
```



# Toolkit Contents
The following can be performed with this suite of tools:  
* highlight UMAP clusters of interest _in situ_ across your FOV of choice
* create publication ready plots
* tally the number of cells in each Xenium Assay Seurat object
* merge and analyze data across multiple slices, via [UMAP dimensionality reduction](https://www.nature.com/articles/nbt.4314) and applying [cluster-based algorithms](https://www.tandfonline.com/doi/full/10.1080/15476286.2020.1728961)
* find marker genes for unique UMAP clusters
* create Box plots for cluster-specific marker genes (using sequencing depth-corrected or raw counts)
* calculate the midline of any group of cells _in situ_
* compute gene expression as a function of distance away from a computed midline
* perform a gene expression gradient analysis (via computing 1-Dimensional UMAP embeddings values for cells)

# Function List
`HighlightCluster`
* highlight cluster(s) of choice in an Xenium object  
{_cluster of interest labelled "this cluster" & all other clusters in the object labelled "all else"_}

`HighlightCells`
* highlight subset cells from a Xenium object of choice  -- useful when subsetting cells for downstream analysis  
{_cells from subset object labelled "these cells" & all other cells labelled "all else"_}

`XeniumBoxPlot`
* plot boxplots of given genes per cluster using sequencing-depth corrected counts  
  
`XeniumBoxPlotRaw`
* plot boxplots of given genes per cluster using raw Xenium counts  

additional functions provided by Mark S. Cembrowski  **update these**  
`getCentre`

`getDistanceToLine`

> [!IMPORTANT]  
> **Prior to running any 1D UMAP plots (i.e. plotUMAP1_inSituMidline) UMAP_1 embeddings must be computed.**

User can compute UMAP_1 embeddings via:  
```
#subset your cluster(s) of interest
AD = xen_atn_subregions %>% subset(idents = c('5', '6'))
#(optional) visualize and check cells included in subset 
Seurat::ImageDimPlot(AD, cols = "polychrome", size = 2, fov = 'X1fov') 

#compute 1D UMAP embeddings
AD<- Seurat::RunPCA(AD, npcs = 30, features = rownames(AD))
#choose number of dims for UMAP accordingly
Seurat::ElbowPlot(AD, ndims = 30, reduction = "pca")
AD <- Seurat::RunUMAP(AD, dims = 1:15, n.components = 1)
```
*note*: these functions are compatible with Seurat, for more information see [https://satijalab.org](https://satijalab.org/seurat/articles/install_v5)

[back to top](#top)
