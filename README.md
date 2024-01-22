# Xenium Spatial Analysis R Toolkit  
<img src="https://github.com/MargoKapustina/Xenium-spatial-tools/assets/129800017/e3f99642-03a7-4b90-ae9e-a85eab71bb97" width="20%"></img>    
This repository contains R tools and workflows to spatially analyze cell-types using single-cell spatial transcriptomic *Xenium* data.   

## Table of Contents
1. [Prep Environment](#prep-environment-installing-seurat)
     -deprecate-
3. [Toolkit Contents](#toolkit-contents)
4. [Function List](#function-list)


# Prep Environment: installing Seurat
These tools are compatible with Seurat V5. For more information, see https://satijalab.org/seurat/articles/install_v5.   
To install and load package:
```
install.packages('Seurat')
library(Seurat)
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
