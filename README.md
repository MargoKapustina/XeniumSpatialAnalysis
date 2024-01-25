# Xenium Spatial Analysis R Toolkit  
<img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/e3f99642-03a7-4b90-ae9e-a85eab71bb97" width="20%"></img>    
This repository contains R tools and workflows to spatially analyze cell-types using single-cell spatial transcriptomic *Xenium* data.   

## Table of Contents
1. [Installation](#installation)
2. [Tutorial](#tutorial)
3. [Toolkit Contents](#toolkit-contents)
4. [Function List](#function-list)


# Installation
To install and load package:
```R
devtools::install_github("MargoKapustina/XeniumSpatialAnalysis")
library(XeniumSpatialAnalysis)
```
When updating to a newer version of the repo:
```R
#remove old version
remove.packages(XeniumSpatialAnalysis)  

#reinstall from here or from the cembrowskilab/RUHi github  
devtools::install_github("MargoKapustina/XeniumSpatialAnalysis")
```

# Tutorial 
 
## 1. Read in Data ##
### Read in raw Xenium data ### 
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
### Merge data across slices and animals ### 
Read in data from `.rds` files, add sample `metadata` to differentiate slices, and `merge` slices. 
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
### 3a. Analyze merged data ### 
Analyze data across multiple slices, via [UMAP dimensionality reduction](https://www.nature.com/articles/nbt.4314) and applying [cluster-based algorithms](https://www.tandfonline.com/doi/full/10.1080/15476286.2020.1728961). Continue until you have desired **UMAP embeddings**. 
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

#find the number of cells in your object
dim(xenexc) #75058 - excitatory neuronal cells
```
To highlight a cluster in your Xenium object, use `HighlightCluster()` and specify:
* `obj`: Xenium object
* `cluster_id`: Identity of cluster to highlight
* `FOV`: FOV(s) of choice to vizualise cells within
> Following are __optional__ to specify:
> * `size`: Size of the highlighted cells in the plot  
> * `alpha_value`: Alpha (transparency) value of the highlighted cells  
> * `color_palette`: Color palette for the plot  
> * `save_plot`: Option to save plot as .pdf in working directory (TRUE, FALSE) 
```R
highlightCluster(obj = xenexc, cluster_id = 1, FOV = c('X1fov', 'X2fov'))
```
### 3a. Subset clusters of choice for downstream analysis ### 
Subset clusters from a Xenium object with `subset()`.
```R
#subset clusters
xen_atn = xenexc %>% subset(idents = c("2", "6"))
```
#### Next, to highlight your subset object within the original object, use `HighlightCells()` and specify: #####
* `highlight_obj`: the object you want to highlight (ie subset cells)
* `within_obj`: the object you want to highlight your object within (ie cells plotted but not highlighted; original object)
* `FOV`: FOV of choice to vizualise cells within  
> __Optional__ to specify:
> * `size`: Size of the highlighted cells in the plot  
> * `alpha_value`: Alpha (transparency) value of the highlighted cells  
> * `color_palette`: Color palette for the plot  
> * `save_plot`: Option to save plot as .pdf in working directory (TRUE, FALSE) 
```R
#highlight subset of cells within original object (or any smaller object within a larger object)
highlightCells(highlight_obj = xenatn, within_obj = xenexc)
```
> [!TIP]
> #### Repeat these steps until you are satisfied with your final subset of cells. ####
(ie example ATN cells!)

## 4. Begin spatial analysis! ##
#### 4a. Create marker gene boxplots for each cluster #### 
Create boxplots for gene of interest across all clusters, using `XeniumBoxPlot()` or `XeniumBoxPlotRaw()`
```R
#create vector of genes
my_genes= c('C1ql2','Slc17a7', 'Gng13')

#generate boxplots using sequencing depth corrected counts (these are SCT$counts)
BoxPlots = XeniumBoxPlot(object = atn, genes = my_genes)

#or using raw counts (these are Xenium$counts)
BoxPlotsRaw = XeniumBoxPlot(object = atn, genes = my_genes)
```

### 4b. Assess spatial gradients in gene expression ###
Create a 1-dimensional UMAP, plot the corresponding histogram of UMAP_1 embedding values and the UMAP_1 embedding values in situ within specified FOV. Additionally, compute the midline intersecting cells within specified FOV.
__note:__ please compute UMAP_1 embeddings beforehand! :
```R
#subset cells in specific clusters (see 3a)
AD = xen_atn_subregions %>% subset(idents = c('5', '6'))
#visualize the subset
ImageDimPlot(AD, cols = "polychrome", size = 2, fov = 'X1fov')
#highlight cells that we subset
highlightCells(highlight_obj = AD, within_obj = xen_atn_subregions)

#compute 1D UMAP embeddings
AD<- RunPCA(AD, npcs = 30, features = rownames(AD))
#choose numbers of dims based on elbowplot
ElbowPlot(AD, ndims = 30, reduction = "pca") 
AD <- RunUMAP(AD, dims = 1:15, n.components = 1)
```
### 4c. Assess spatial gradients in gene expression ###
Run `plotUMAP1_inSituMidline()`, making sure to check how your spatial midline looks. If it looks off - please adjsut `degs` parameter. Specify:
* `object` a Xenium object
* `gene`s Character vector of gene names to fetch expression data for
* `FOV` FOV to extract coordinates from
* `degs` User-defined angle, defined in degrees, that intersects centrepoint in midline computation
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)
```R
#Plot 1D UMAP, corresponding histogram of UMAP_1 embedding values, UMAP_1 embedding values in situ (within specified FOV), and compute the midline intersecting cells (within specified FOV).
#note: please compute UMAP_1 embeddings beforehand.
AD_df_midline = plotUMAP1_inSituMidline(object = AD, FOV = 'fov', degs = 65, save_plot = TRUE)
```

You can also create the above plots without computing the spatial midline using `plotUMAP1_inSitu()` Specify:
* `object` a Xenium object
* `FOV` FOV to plot UMAP_1 embedding values in
> * `EmbeddingsPlotTitle` title for UMAP_1 Embeddings plot
> * `HistogramPlotTitle` title for UMAP_1 Embeddings Histogram plot
> * `inSituPlotTitle `title for UMAP_1 Embeddings in situ plot
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)
```R
#Plot 1D UMAP, corresponding histogram of UMAP_1 embedding values, UMAP_1 embedding values in situ (within specified FOV)
#note: please compute UMAP_1 embeddings beforehand.
AD_df = plotUMAP1_inSitu(object = AD, FOV = 'fov', save_plot = TRUE)
```

### 5a. Assess spatial gradients in expression of individual genes ###
First generate Gene Expression vs Midline data for indivudal FOVs using `getExpressionvsMidline()`. When using this function, please specify:
* `object` a Xenium object
* `genes` Character vector of gene names to fetch expression data for
* `FOV` FOV to extract coordinates from
* `degs` User-defined angle, defined in degrees, that intersects centrepoint in midline computation
```R
#generate gene expression vs spatial midline data for individual FOVs, for any gene(s) in your assay
geneExpression_X1fov <- getExpressionvsMidline(AD, FOV = "X1fov", genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'), degs = 45)
geneExpression_fov <- getExpressionvsMidline(AD, FOV = "fov", genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'), degs = 50)
```
Then, use `plotGeneExpressionVsMidline()` to create a pooled dataframe containing Cell IDs, coordinates (X,Y), gene expression counts for specififed gene(s), and computed distance away from spatial midline of each cell, across all FOVs in data. Running `getExpressionvsMidline()` will also plot gene expression data for cells and their distance away from Spatial Midline across _multiple_ FOVs. When using this function, please specify:
* `geneExpressionData` List of dataframes generated with getExpressionvsMidline() across multiple FOVs
* `genes` Character vector of gene names to fetch expression data for
> __Optional__ to specify:
> * `binNumber` Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth (default: 7)
> * `binwidth` Width of bins for gene expression averaging (lines) Suggested binwidth = 40 for 40micron bins (default: 40)
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)
```R
#pool data and plot gene expression vs spatial midline data for multiple genes, across multiple FOVs
#lines are gene expression values (averaged gene counts per binned cells) & histogram is number of cells in each bin 
pooled_GeneExpression = plotGeneExpressionVsMidline(geneExpressionData = list(geneExpression_fov, geneExpression_fX1ov), genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'))
```
> You can also run the function on just a single genExpression dataframe: 
```R
#pool data and plot gene expression vs spatial midline data for multiple genes, for one FOV
#lines are gene expression values (averaged gene counts per binned cells) & histogram is number of cells in each bin 
pooled_GeneExpression = plotGeneExpressionVsMidline(geneExpressionData = list(geneExpression_fov), genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'))
```
### Selecting `genes`: ###  
> It's recommended that you specifiy the same `genes` for `plotGeneExpressionVsMidline()` as `getExpressionvsMidline()`. If you would like to get gene expressiond data for more genes, run `getExpressionvsMidline()` with your new desired genes, and then run `plotGeneExpressionVsMidline()` after. If you keep the dataframe name the same, running the function will overwrite the previous gene expression dataframe stored in your environment.

In a similar fashion to plotting gene expression data as a function of distance away from spatial midline, you can also plot the UMAP_1 embeddings as a function of distance away from midline. To do this, use `getUMAP1_MidlineData()` then `plotUMAP1_vsMidline()`.
When using `getUMAP1_MidlineData()` please specify:
* `object` a Xenium object
* `FOV` FOV to extract coordinates and UMAP_1 embeddings values from
> __Optional__ to specify:
> * `binNumber` Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth (default: 7)
> * `binwidth` Width of bins for gene expression averaging (lines) Suggested binwidth = 40 for 40micron bins (default: 40)
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)  

When using `plotUMAP1_vsMidline()` please specify:
* `UMAP1MidlineData` List of dataframes generated with getUMAP1_MidlineData() across multiple FOVs
> __Optional__ to specify:
> * `binNumber` Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth (default: 7)
> * `binwidth` Width of bins for gene expression averaging (lines) Suggested binwidth = 40 for 40micron bins (default: 40)
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)  
```R
#generate UMAP_1 embedding values vs spatial midline data for individual FOVs, specifying degrees to define spatial midline for each individual FOV
#will also plot coloured coordinates according to UMAP_1 embeddings and histogram of UMAP_1 embeddings values per binned distance for single FOV specified
##note: please compute UMAP_1 embeddings beforehand
UMAP1_midline_data_fov <- getUMAP1_MidlineData(AD, FOV= c('fov'), degs = 65)
UMAP1_midline_data_X1fov <- getUMAP1_MidlineData(AD, FOV= c('X1fov'), degs = 45)

#pool data and plot histogram of UMAP_1 embeddings values per binned distance across multiple FOVs
pooled_UMAP1_vsMidline = plotUMAP1_vsMidline(list(UMAP1_midline_data_fov, UMAP1_midline_data_X1fov))
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
