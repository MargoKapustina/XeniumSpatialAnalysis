# Xenium Spatial Analysis Toolkit  
🆇🅴🅽🅸🆄🅼 🆂🅿🅰🆃🅸🅰🅻 🅰🅽🅰🅻🆈🆂🅸🆂   
This repository contains tools and workflows to spatially analyze cell-types using single-cell spatial transcriptomic Xenium data.   

## Table of Contents
1. [Prep Environment](#prep-environment-installing-seurat)
2. [Toolkit Contents](#toolkit-contents)
3. [Function List](#function-list)
   
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


note: prior to running any UMAP_1 embedding plots (plotUMAP1_inSituMidline) UMAP_1 embeddings needs to be computed.
User can run UMAP embeddings computation first via:
```
#subset your cluster(s) of interest
AD = xen_atn_subregions %>% subset(idents = c('5', '6')) #gets cells in these clusters
#visualize and check subset
ImageDimPlot(AD, cols = "polychrome", size = 2, fov = 'X1fov') 

#compute 1D UMAP embeddings
AD<- RunPCA(AD, npcs = 30, features = rownames(AD))
#choose number of dims for UMAP accordingly
ElbowPlot(AD, ndims = 30, reduction = "pca")
AD <- RunUMAP(AD, dims = 1:15, n.components = 1)
```
[back to top](#top)

**want a short description of the scripts that i haave. Check the folder.

allows anyone to create simple boilerplate documentation for any bioinformatics software, using guidelines that have been co-created by the Australian BioCommons and members of the community.

This template repository contains a set of guidelines for documenting bioinformatics software (e.g. tools and workflows).
The initial guideline version uploaded to GitHub was informed by current documentation practices and structures used in the GitHub community.
Subsequent versions have been modified and updated with input from Australian BioCommons engagements with infrastructure partners and the bioinformatics community.
Typical files are included, such as a LICENSE, CITATION.cff and change_log.md
These guidelines will be further developed as needed to meet the requirements of the Australian BioCommons community.

