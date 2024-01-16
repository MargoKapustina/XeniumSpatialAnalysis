# Xenium Spatial Analysis Toolkit
This repository contains tools and workflows to spatially analyze cell-types using Xenium data.   
🆇🅴🅽🅸🆄🅼 🆂🅿🅰🆃🅸🅰🅻 🅰🅽🅰🅻🆈🆂🅸🆂  
## Table of Contents
1. [Prep Environment](#prep-environment-installing-seurat)
2. [Toolkit Contents](#toolkit-contents)
   
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
* Create publication ready plots
* tally the number of cells in each Xenium Assay Seurat object
* merge and analyze data across multiple slices, via [UMAP dimensionality reduction](https://www.nature.com/articles/nbt.4314) and applying [cluster-based algorithms](https://www.tandfonline.com/doi/full/10.1080/15476286.2020.1728961)
* find marker genes for unique UMAP clusters
* create Box plots for cluster-specific marker genes (using sequencing depth-corrected or raw counts)
* calculate the midline of any group of cells _in situ_
* compute gene expression as a function of distance away from a set midline
* perform a gene expression gradient analysis (via computing 1-Dimensional UMAP embeddings values for cells)

**want a short description of the scripts that i haave. Check the folder.

allows anyone to create simple boilerplate documentation for any bioinformatics software, using guidelines that have been co-created by the Australian BioCommons and members of the community.

This template repository contains a set of guidelines for documenting bioinformatics software (e.g. tools and workflows).
The initial guideline version uploaded to GitHub was informed by current documentation practices and structures used in the GitHub community.
Subsequent versions have been modified and updated with input from Australian BioCommons engagements with infrastructure partners and the bioinformatics community.
Typical files are included, such as a LICENSE, CITATION.cff and change_log.md
These guidelines will be further developed as needed to meet the requirements of the Australian BioCommons community.

