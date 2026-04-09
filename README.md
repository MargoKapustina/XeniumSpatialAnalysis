# Xenium Spatial Analysis: R Toolkit  
<p align="center">
 <img src="https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/d776c089-0f7e-4400-9d27-b6ee8dce4125" width="62%"> 
</p>

<p align="center">
  <img src="https://img.shields.io/badge/stable-1.0.1-brightgreen" />
  <img src="https://img.shields.io/badge/dev-1.1.0-blue" />
</p>

This repository contains R tools and workflows for spatial analysis of cell-types using single-cell spatial transcriptomic Xenium data. To know more about Xenium technology, check out the [Xenium workflow](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest).    
***
We will use the anterodorsal nucleus of the anterior thalamic nuclei (ATN) within two brain slices from [our Cell Reports paper](https://www.cell.com/cell-reports/fulltext/S2211-1247(24)00170-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124724001700%3Fshowall%3Dtrue) to demonstrate how to examine gene expression gradients. You can also check out the raw data and full Seurat processed object at our [GEO portal](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255953), or download just the anterodorsal nucleus of the ATN in two slices [here](https://osf.io/zae5v/?view_only=27b80485b12441478755fe61613743e3) in `.rds` file format to follow the tutorial. 

_April 6th 2026: New functions (interactiveCellSelector and plotMetadataFrequency) have been added along with tutorials!_
* `interactiveCellSelector` allows you to plot out Xenium object cell spatial coordinates, interactively lasso (i.e., drag-to-select) a region of interest, and save selected cell IDs in your environment.
* `plotMetadataFrequency` allows you to generate a cluster frequency or proportion plot stratified by metadata directly from a Seurat object.  

_July 9th 2024: New functions (analyzeLayer and object_FOV_to_coordinates) have been added along with tutorials!_
* `object_FOV_to_coordinates` allows you to extract cell ID, cluster ID, cell coordinates, etc. for cells from a specified FOV in your Xenium object.
* `analyzeLayer` allows you to compute cell normalized distances along along a cortical layer, and away from a boundary that you define to examine superficial-deep differences in cell-type composition or spatial organization.

## Overview
1. [Toolkit Contents](#toolkit-contents)
2. [Citation](#citation)
3. [Installation](#installation)
4. [Tutorial](#tutorial)
5. [Function List](#function-list)

# Toolkit Contents
The following can be performed with this suite of tools:  
* create publication ready plots
* merge and analyze data across multiple slices, via [UMAP dimensionality reduction](https://www.nature.com/articles/nbt.4314) and applying [cluster-based algorithms](https://www.tandfonline.com/doi/full/10.1080/15476286.2020.1728961)
* highlight UMAP clusters of interest _in situ_ across selected FOV(s)
* tally the number of cells in each Xenium Assay Seurat object
* find marker genes for unique UMAP clusters
* create boxplots for cluster-specific marker genes (using sequencing depth-corrected or raw counts)
* calculate the midline of any group of cells _in situ_
* compute gene expression as a function of distance away from a computed midline
* perform a gene expression gradient analysis (via computing 1-Dimensional UMAP embeddings values for cells)
* compute UMAP_1 embeddings as a function of distance away from a computed midline
* spatially subset cell IDs from a Xenium object with an interactive lasso tool
* ... and more!

# Citation
If you find this Xenium Toolkit useful, please cite our paper:
```text
Kapustina, M., Zhang, A.A., Tsai, J.Y.J., Bristow, B.N., Kraus, L., Sullivan, K.E., Erwin, S.R., Wang, L., Stach, T.R., Clements, J., et al. (2024). The cell-type-specific spatial organization of the anterior thalamic nuclei of the mouse brain. Cell Reports 43, 113842–113842. https://doi.org/10.1016/j.celrep.2024.113842.
```

# Installation 
> **Please ensure that your R version >= 4.3.2.**     

To install and load package:
```R
devtools::install_github("MargoKapustina/XeniumSpatialAnalysis")
library(XeniumSpatialAnalysis)
```
When updating to a newer version of the repo:
```R
#remove old version
remove.packages("XeniumSpatialAnalysis")  

#reinstall from here or from the cembrowskilab/XeniumSpatialAnalysis github  
devtools::install_github("MargoKapustina/XeniumSpatialAnalysis")
```

# Tutorial 
 
## 1. Read in Data ##
### Read in raw Xenium data ### 
Read in raw Xenium data using Seurat's `LoadXenium()` function. Make sure to set unique FOV names for each slice for easier downstream processing. The code here also removes cells with zero counts, and saves data as .rds file in a specified directory. 

For each replicate assign a unique FOV for easier downstream processing:
```R
#Author: Margarita Kapustina, UBC, 2023
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
#load packages
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

#read in data from saved .rds files and merge all slices  
data.dir = "myobjpath1.rds"
S1_ant = readRDS(data.dir)
data.dir = "myobjpath2.rds"
S1_int = readRDS(data.dir)
data.dir = "myobjpath3.rds"
S1_post = readRDS(data.dir)

#add sample metadata
S1_ant@meta.data$sample = 'S1_ant'
S1_int@meta.data$sample = 'S1_int'
S1_post@meta.data$sample = 'S1_post'

#merge slices
merge = merge(S1_ant, S1_int)
merge = merge(merge, S1_post)
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
#visualize subset of cells in specified FOV
ImageDimPlot(xenexc, fov = 'X1fov')

#run analysis
xenexc <- RunPCA(xenexc, npcs = 30, features = rownames(xenexc))
#select your number of dims based on elbow plot
     #fyi elbow plot: a ranking of principle components based on the
     #percentage of variance explained by each one (see: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
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
highlightCluster(obj = xenexc, cluster_id = c('2', '6'), size = 2, FOV = c('X1fov', 'fov'))
```
<p align="center">
  <img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/303f28c5-da7b-47de-acd7-9e8113497a22">
</p>   

### 3b. Subset clusters of choice for downstream analysis ### 
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
highlightCells(highlight_obj = xen_atn, within_obj = xenexc, size = 1)
```
<p align="center">
  <img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/89a2bc2e-0975-4b93-88b9-a0e7d6b5cdb3"
</p>   

#### Repeat these steps until you are satisfied with your resulting subset of cells. ####

## 4. Prep your spatial gradient analysis ##
#### 4a. Create marker gene boxplots for each cluster #### 
Create boxplots for gene of interest across all clusters, using `XeniumBoxPlot()` or `XeniumBoxPlotRaw()`
```R
#create vector of genes
my_genes= c('C1ql2','Foxp1', 'Gng13')

#generate boxplots using sequencing depth corrected counts (these are SCT$counts)
BoxPlots = XeniumBoxPlot(object = xen_atn, genes = my_genes)
#or using raw counts (these are Xenium$counts)
BoxPlotsRaw = XeniumBoxPlotRaw(object = xen_atn, genes = my_genes)
```
<p align="center">
  <img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/b4252ee9-ceb6-4908-90bc-59fc7075fbd1"
</p>   

### 4b. Visualize clusters chosen for analysis in space ###
Before performing your spatial gradient analysis or compute the spatial midline, visualize your cells in space with `highlightCells()`.
```R
#subset cells in specific clusters (see 3a)
AD = xen_atn_subregions %>% subset(idents = c('5', '6'))
#visualize the subset
ImageDimPlot(AD, cols = "polychrome", size = 2, fov = 'X1fov')
#highlight cells that we subset
highlightCells(highlight_obj = AD, within_obj = xen_atn_subregions)
```
<p align="center">
<img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/cd23ffa2-8034-4d6e-9bb6-8f842c6f2fe3"
</p>   

> [!IMPORTANT]  
> **Prior to running any 1D UMAP plots (i.e. plotUMAP1_inSituMidline) UMAP_1 embeddings must be computed.**  
> <details>
> <summary>User can compute UMAP_1 embeddings using Seurat </summary>
> <pre>
> #subset your cluster(s) of interest
> AD = xen_atn_subregions %>% subset(idents = c('5', '6'))
>
>  #(optional) visualize and check cells included in subset 
> Seurat::ImageDimPlot(AD, cols = "polychrome", size = 2, fov = 'X1fov') 
>
>  #load required package
> library(Seurat)
>  #compute 1D UMAP embeddings
> AD<- Seurat::RunPCA(AD, npcs = 30, features = rownames(AD))
>  #choose number of dims for UMAP accordingly
> Seurat::ElbowPlot(AD, ndims = 30, reduction = "pca")
> AD <- Seurat::RunUMAP(AD, dims = 1:15, n.components = 1)
>
> For more more info on Seurat functions see (https://satijalab.org/seurat/articles/install_v5) </pre>
> </details>

## 5. Perform your spatial gradient analysis ##
### 5a. Perform spatial gene expression analysis using UMAP_1 embeddings ###
Here, we will:    
i. Create a 1-dimensional UMAP   
ii. Plot the corresponding histogram of UMAP_1 embedding values & the UMAP_1 embedding values *in situ* within specified FOV   
iii. Compute the midline intersecting cells within specified FOV     

  
First, run `plotUMAP1_inSituMidline()`, making sure to check how your spatial midline looks. If it looks off - please adjsut `degs` parameter. When using the function specify:
* `object` a Xenium object
* `FOV` FOV to extract coordinates from
* `degs` User-defined angle, defined in degrees, that intersects centrepoint in midline computation
> * `EmbeddingsPlotTitle` title for UMAP_1 Embeddings plot
> * `HistogramPlotTitle` title for UMAP_1 Embeddings Histogram plot
> * `inSituPlotTitle `title for UMAP_1 Embeddings in situ plot
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)
```R
#Plot 1D UMAP, corresponding histogram of UMAP_1 embedding values, UMAP_1 embedding values in situ (within specified FOV), and compute the midline intersecting cells (within specified FOV).
#note: please compute UMAP_1 embeddings beforehand.
AD_df_midline = plotUMAP1_inSituMidline(object = AD, FOV = 'X1fov', degs = 45, save_plot = TRUE)
```
<p align="center">
<img width = '75%'src="https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/b0ba20b7-e03d-4bb2-9b23-e762aef10ed1"
</p>   


You can also create the above plots without computing the spatial midline using `plotUMAP1_inSitu()` When using the function specify:
* `object` a Xenium object
* `FOV` FOV to plot UMAP_1 embedding values in
> * `EmbeddingsPlotTitle` title for UMAP_1 Embeddings plot
> * `HistogramPlotTitle` title for UMAP_1 Embeddings Histogram plot
> * `inSituPlotTitle `title for UMAP_1 Embeddings in situ plot
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)
```R
#Plot 1D UMAP, corresponding histogram of UMAP_1 embedding values, UMAP_1 embedding values in situ (within specified FOV)
#note: please compute UMAP_1 embeddings beforehand.
AD_df = plotUMAP1_inSitu(object = AD, FOV = 'X1fov', save_plot = TRUE)
```
<p align="center">
<img width = '77%' src="https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/6bc805ec-c823-41f2-9542-07ff80c52184"
</p>   


### 5b. Perform spatial gene expression analysis using individual genes ###
First generate Gene Expression vs Midline data for individual FOVs using `getExpressionvsMidline()`. Make sure to visually assess whether your midline intersects your cells appropriately. When using this function, please specify:
* `object` a Xenium object
* `genes` Character vector of gene names to fetch expression data for
* `FOV` FOV to extract coordinates from
* `degs` User-defined angle, defined in degrees, that intersects centrepoint in midline computation
```R
#generate gene expression vs spatial midline data for individual FOVs, for any gene(s) in your assay
geneExpression_X1fov <- getExpressionvsMidline(AD, FOV = "X1fov", genes = c("Epha1", "Cntnap4", 'Pcp4', 'Mdga1'), degs = 45)
geneExpression_fov <- getExpressionvsMidline(AD, FOV = "fov", genes = c("Epha1", "Cntnap4", 'Pcp4', 'Mdga1'), degs = 60)
```

<p align="center">
<img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/d9741824-cd75-40af-82af-9d2beeb74410" width="40%"></img>
<img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/8b90de9d-2c71-4efe-9a17-65cf9f6b4af2" width="27%"></img>   
</p>

***
Then, use `plotGeneExpressionVsMidline()` to create a pooled dataframe containing Cell IDs, coordinates (X,Y), gene expression counts for specififed gene(s), and computed distance away from spatial midline of each cell, across all FOVs in data. Running `getExpressionvsMidline()` will also plot gene expression data for cells and their distance away from Spatial Midline across _multiple_ FOVs. When using this function, please specify:
* `geneExpressionData` List of dataframes generated with getExpressionvsMidline() across multiple FOVs
* `genes` Character vector of gene names to fetch expression data for
> __Optional__ to specify:
> * `binNumber` Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth (default: 7)
> * `binwidth` Width of bins for gene expression averaging (lines) Suggested binwidth = 40 for 40micron bins (default: 40)
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)
> * `xlim` Define x-axis limits as vector (default: -500 to 500)
```R
#pool data and plot gene expression vs spatial midline data for multiple genes, across multiple FOVs
#lines are gene expression values (averaged gene counts per binned cells) & histogram is number of cells in each bin 
pooled_GeneExpression = plotGeneExpressionVsMidline(geneExpressionData = list(geneExpression_fov, geneExpression_X1fov), genes = c("Epha1", "Cntnap4", 'Pcp4', 'Mdga1'))
```
<p align="center">
<img width = '77%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/6d1097da-6be0-4a7c-8d41-b43f36e0ac55"
</p>     


<details>
<summary>You will also get a readout of the colour legend per genes and the number of cells per bin when running this function</summary>
<pre>Colour Scheme (Gene : colour):
Epha1 : red 
Cntnap4 : green 
Pcp4 : pink 
Mdga1 : black 
[1] "Generating plot..."
Number of cells in each bin (note: plotted are counts/100):
Beginning on left-most bin...
Bin: 1 , Count: 11 
Bin: 2 , Count: 28 
Bin: 3 , Count: 66 
Bin: 4 , Count: 78 
Bin: 5 , Count: 29 
Bin: 6 , Count: 21 
Bin: 7 , Count: 1 </pre>
</details>  

> You can also run the function on just a single genExpression dataframe: 
```R
#pool data and plot gene expression vs spatial midline data for multiple genes, for one FOV
#lines are gene expression values (averaged gene counts per binned cells) & histogram is number of cells in each bin 
pooled_GeneExpression = plotGeneExpressionVsMidline(geneExpressionData = list(geneExpression_fov), genes = c("Epha1", "Cntnap4", 'Pcp4', 'Mdga1'))
```
#### _note on selecting `genes`_: ####
> It's recommended that you specify the same `genes` for `plotGeneExpressionVsMidline()` as `getExpressionvsMidline()`. If you would like to get gene expression data for more genes, run `getExpressionvsMidline()` with your new desired genes, and then run `plotGeneExpressionVsMidline()` after. If you keep the dataframe name the same, running the function will overwrite the previous gene expression dataframe stored in your environment.
***
In a similar fashion to plotting gene expression data as a function of distance away from spatial midline, you can also plot the UMAP_1 embeddings as a function of distance away from midline. To do this, use `getUMAP1_MidlineData()` then `plotUMAP1_vsMidline()`.   
When using `getUMAP1_MidlineData()` please specify:
* `object` a Xenium object
* `FOV` FOV to extract coordinates and UMAP_1 embeddings values from
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
```
<p align="center">
<img src = "https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/fce46bba-5e8d-4f69-b0a1-1f69bd24e390"
</p>   

When using `plotUMAP1_vsMidline()` please specify:
* `UMAP1MidlineData` List of dataframes generated with getUMAP1_MidlineData() across multiple FOVs
> __Optional__ to specify:
> * `binNumber` Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth (default: 7)
> * `binwidth` Width of bins for gene expression averaging (lines) Suggested binwidth = 40 for 40micron bins (default: 40)
> * `save_plot` Option to save plot as .eps in working directory (TRUE, FALSE)  
 
```R
#pool data and plot histogram of UMAP_1 embeddings values per binned distance across multiple FOVs
pooled_UMAP1_vsMidline = plotUMAP1_vsMidline(list(UMAP1_midline_data_fov, UMAP1_midline_data_X1fov))
```
<p align="center">
<img src = "https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/2643f1a1-d079-4f9c-bde7-5888cdb293e5"
</p>

### New Functions:
To examine whether your cells exhibit superficial-deep transcriptomic differences, the following examples will guide you through extracting cell coordinates with relevant information and computing distance along a cortical layer and cell distances away from a boundary that you define.

First, use `object_FOV_to_coordinates()` to create a dataframe with cell ID coordinates and cluster IDs from a Xenium Seurat object. This function will also allow you to define a boundary for each FOV, that you specifiy if you are interested in assessing whether cell cluster IDs differ above, and below, a regional boundary that you define. The function will return a dataframe with cell ID, coordinates, cluster ID, and a plot with UMAP_1 embedding values and the boundary that you define. Please make sure that you have already run dimensionality reduction and clustering on your Xenium object before running this function.

When using `object_FOV_to_coordinates()` please specify:
* `object` Xenium object to extract cell coordinates from
* `thisFOV` Name of the FOV to extract coordinates from
* `threshold_y` Boundary that you define to compute cells that are above and below this
* `angle_adjust` Option to adjust the angle of your cell coordinates (TRUE, FALSE)
* `theta_deg` Specify the degrees you wish to rotate your slice by
* `flip_x_coordinates` Option to flip your slice in the horizontal plane (TRUE, FALSE)

```R
distanceData_X11fov = object_FOV_to_coordinates(object, thisFOV = 'X11fov', threshold_y = 1911, 
                                                angle_adjust = TRUE, 
                                                theta_deg = 20, 
                                                flip_x_coordinates = TRUE)
```

<p align="center">
<img width = '60%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/5c6e6171-f5ef-460c-a28f-dfd827ce05a2"
</p>

Once you have run object_FOV_to_coordinates() on all slices, merge your FOV-specific dataframes using rbind(), and compute the normalized distance away from the boundary for each of the slices.

```R
#merge all distancedistanceData
distance_data_all = rbind(distanceData_X2fov, distanceData_X3fov)
distance_data_all = rbind(distance_data_all, distanceData_X4fov)
distance_data_all = rbind(distance_data_all, distanceData_X11fov)

#add normalized distance
distance_data_all$normalized_distance = (distance_data_all$y - distance_data_all$threshold)
```
Now, use `analyzeLayer()` to compute the cell normalized distances along the cortical layer and away from a boundary that you define. To define your boundary, click sequential points in your plot screen, and hit `ESC` when complete. The boundary defined will be showed in the plot screen.
When using `analyzeLayer()` please specify:
* `thisFOV` Name of our FOV to extract coordinates from

__Make sure that you have your merged distanceData stored as distance_data_all__
```R
#make sure that you have merged distanceData stored as distance_data_all
#if you only have one FOV, rename the dataframe to distance_data_all

#save your working directory to where plot output will be saved
setwd("/myPlots")
cellDataX11fov <- analyzeLayer("X11fov")
dev.copy2eps(file = 'X11fovAnalysis.eps')
```
<p align="center">
<img width = '60%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/1ef250f0-bbdc-47b1-ab8c-3d7e00d2d47c"
</p>

_Upon running this function, several plots will be saved. These plots show:_
> * normalized distance along the cortical layer
> * normalized distance to boundary of the cortical layer
> * histogram of distances along the cortical layer
> * histogram of distances to the boundary of the cortical layer

<p align="center">
<img width = '41%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/d9ef686e-c52f-4340-a322-2fdd794fdd55"></img>
<img width = '40%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/f2d91f44-7ff9-47b0-91c1-4c09a6f0ed64"></img>
<img width = '40%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/2001adea-e51f-43c7-87fb-a68937b2d861"></img>
<img width = '40%' src = "https://github.com/MargoKapustina/XeniumSpatialAnalysis/assets/129800017/aea3dca3-ae7f-4b39-b4b8-eb00e3ca82d0"></img>
</p>

Next, you can merge FOV-specific cell data into a merged dataframe to compare cell cluster ID and normalized distance along the neocortical layer, or away from the defined boundary
```R
#merge all cellData
allCellData = rbind(cellDataX2fov, cellDataX3fov)
allCellData = rbind(allCellData, cellDataX4fov)
allCellData = rbind(allCellData, cellDataX11fov)
```
--- 
Next, to spatially select cells and retrieve their cell ID, for the purpose of adding metadata or subsetting these cells, the following example will walk you through how to use an interactive lasso function and example use-cases. 

First, before running this function, use `object_FOV_to_coordinates()` to create a dataframe with cell ID coordinates and cluster IDs from a Xenium Seurat object.

Use `interactiveCellSelector()` to interactively select cell IDs from a spatial region of interest. This function will launch a Shiny app allowing you to spatially lasso cells. The app displays an interactive spatial plot of cells colored by cluster identity. Users can drag to select cells of interest. Upon selection, click the **"Save Selected Cell IDs"** button, which will save two objects in the global environment: 
* (1) `cellid_<your_fov_label>` : a vector of selected cell IDs (e.g., distanceData_X1fov). This vector can be used to subset cells via: `subsetObj <- subset(myObj, cells = cellid_yourFovLabel)`
* (2) `<your_fov_label>` : a filtered dataframe containing only the selected cells. The resulting object will be save using your specified `fov_label` as its name (e.g., x1fov). 

When using `interactiveCellSelector()` please specify:
* `fov_label` A character string specifying the field of view (FOV) label (e.g., 'x1fov')
* `cell_coords_df` A dataframe containing cell coordinate data, generated using `object_FOV_to_coordinates()`(e.g., 'distanceData_X2fov')
> __Optional__ to specify:
> * `plot_height`: Character string specifying plot height. This controls the vertical display size of the interactive visualization. Must include px units (e.g. "750px").
> * `pt.size`: Numeric value controlling point size in the scatter plot (default 0.8)

```R
#launch interactive cell ID lasso selection
interactiveCellSelector("x2fov", distanceData_X2fov)
# step 1: Lasso cells of interest.
# step 2: click "Save Selected Cell IDs" button
# step 3: close Shiny app window
# Result: cellid_x2fov (selected cell ID vector) and x2fov (filtered selected cell coordinates dataframe) will now be present in your global environment.
# Note: if cell ID vector is missing but filtered coordinate dataframe is present, run cellID_x2fov = x2fov$cell to store cell ID vector manually
```
<p align="center">
 <img width = '80%' src="https://github.com/user-attachments/assets/6d17f339-abf5-4541-83f0-9b88b2f77527" alt="lassoFastExample" />
</p> 


You can use the generated cell ID vector and filtered coordinate dataframe for downstream analysis. Some examples are:
```R
#load seurat if needed
library(Seurat)
## To subset interactiveCellSelector() selected cell IDs via from a Seurat object:
subsetObj <- subset_opt(myObj, cells = cellID_x2fov)

# OPTIONAL: To subset interactiveCellSelector() cell IDs which are not present across all sections (i.e. FOVs)
# first store the modified subset_opt function from: https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat.R  
# then run:
subsetObj <- subset_opt(myObj, cells = cellID_x2fov)

# OPTIONAL: To combine cell ID vector saved via interactiveCellSelector() then subset from a Seurat object:
combinedCellids <- c(distanceData_X1fov, distanceData_X2fov)
subsetObj <- subset(object, cells = combinedCellids)

## To add metadata to cells selected via interactiveCellSelector():
# 1) add variable name of interest to generated df
x2fov$metadataVariable = 'metadataVariable'
# 2) create metadata dataframe using selected cell IDs
metadata <- data.frame(
  metadataVariable = x2fov$metadataVariable, 
  row.names = x2fov$cell)
# 3) add metadata to Seurat object
myObj <- AddMetaData(myObj, metadata = metadata, col.name = "metadataVariable")
```

To generate cluster frequency or proportion bar plots stratified by metadata directly from a Seurat object, you can use `plotMetadataFrequency`. Two plots will be returned. (1) Bar plot with raw cell counts per cluster. (2) A normalized bar plot showing the relative proportion of each metadata group within each cluster.
When using `plotMetadataFrequency` please specify:
* `obj` : your Seurat object
* `cluster_var`: Character. Metadata column containing cluster labels.
* `group_var` Character. Metadata column used for grouping (fill variable).
> __Optional__ to specify:
> * `facet_var` Character. Metadata column to facet plots by. Default is NULL.
> * `colors` Character vector of colors. Used if number of groups is <= length(colors).
```R
plots <- plotMetadataFrequency(obj = myObj, 
    cluster_var = "seurat_clusters", 
    group_var = "sample")

#visualize plots
plots$raw+plots$norm
```
<p align="center">
<img width="893" height="583" alt="PlotFrqMeta" src="https://github.com/user-attachments/assets/cf7d855b-4c8a-4105-9a34-e6fc9bfc943b" />
</p>

```R 
plots <- plotMetadataFrequency(obj = myObj, 
    cluster_var = "seurat_clusters", 
    group_var = "sample",
    facet_var = "age")

#visualize plots
plots$raw+plots$norm
```
<p align="center">
<img width="887" height="586" alt="PlotFrqMeta_facet" src="https://github.com/user-attachments/assets/2c7419d2-d174-4aab-bd92-35f0aa8ad3ff" />
</p>

# Function List
`HighlightCluster`
* highlight cluster(s) of choice in a Xenium object  
{_cluster of interest labelled "this cluster" & all other clusters in the object labelled "all else"_}

`HighlightCells`
* highlight subset cells from a Xenium object of choice  -- useful when subsetting cells for downstream analysis  
{_cells from subset object labelled "these cells" & all other cells labelled "all else"_}

`XeniumBoxPlot`
* plot boxplots of given genes per cluster using sequencing-depth corrected counts  
  
`XeniumBoxPlotRaw`
* plot boxplots of given genes per cluster using raw Xenium counts  

`getExpressionvsMidline`
* generates gene expression data away from the computed spatial midline. Supports one FOV at a time

`plotGeneExpressionvsMidline`
* plots gene expression data for cells and their distance away from Spatial Midline across multiple FOVs
  
`getUMAP1_MidlineData`
* fetches UMAP_1 embedding values for cells, and their distance away from the computed spatial midline. Supports one FOV at a time

`plotUMAP1_vsMidline`
* plots UMAP_1 embedding values for cells and their distance away from Spatial Midline across multiple FOVs

`plotUMAP1inSitu`
* plots 1-dimensional UMAP, the corresponding histogram of UMAP_1 embedding values, and the UMAP_1 embedding values in situ within specified FOV
  
`plotUMAP1inSitu_Midline`
* plots gene expression data for cells and their distance away from Spatial Midline across multiple FOVs

#### New Functions
`object_FOV_to_coordinates`
* plots gene expression data for cells and their distance away from Spatial Midline across multiple FOVs

`analyzeLayer`
* plots gene expression data for cells and their distance away from Spatial Midline across multiple FOVs

`interactiveCellSelector`
* interactive lasso-select Shiny app for spatial selection of cell IDs from regions of interest

`plotMetadataFrequency`
* plots cluster frequencies and proportions stratified by metadata, directly from a Seurat object
   
> Interested in recreating the tutorial from our processed Xenium seurat object?
```R
#read in seurat obj from .rds file (access via GEO portal)
xen_atn_analysis = readRDS('xen_atn_analysis.rds')
#increase resolution "high/fine resolution"
xen_atn_subregions = xen_atn_analysis
xen_atn_subregions <- FindClusters(xen_atn_subregions, resolution = 0.7)
AD = xen_atn_subregions %>% subset(idents = c('5', '6'))

## or read in AD data directly (fovs included: fov, X1fov) access via OSF portal
AD = readRDS('exampleAD.rds')
```



[back to top](#top)
