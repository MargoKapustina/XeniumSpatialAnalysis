# Xenium Spatial Analysis Toolkit
🆇🅴🅽🅸🆄🅼 🆂🅿🅰🆃🅸🅰🅻 🅰🅽🅰🅻🆈🆂🅸🆂
This repository contains tools and workflows to spatially analyze cell-types using Xenium data.

[slice-merged1fov.pdf](https://github.com/MargoKapustina/Xenium-tools/files/13956260/slice-merged1fov.pdf)

# Installing Seurat
These tools are compatible with Seurat V5. For more information, see https://satijalab.org/seurat/articles/install_v5. 
To install and load package, use:
```
install.packages('Seurat')
library(Seurat)
```
# Toolkit!
The following can be performed with this suite of tools:
__general__
* highlight UMAP clusters of interest _in situ_ across your FOV of choice
* Create publication-ready plots
* tally the number of cells in each Xenium Assay Seurat object
* merge and analyze data across multiple slices, via UMAP dimensionality reduction and applying cluster-based algorithms
* find marker genes for unique UMAP clusters
* create Box plots for cluster-specific marker genes (using depthc-corrected or raw counts)
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

