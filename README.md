# **SpatialPCA** 

SpatialPCA is a spatially aware dimension reduction method that aims to infer a low dimensional representation of the gene expression data in spatial transcriptomics. SpatialPCA builds upon the probabilistic version of PCA, incorporates localization information as additional input, and uses a kernel matrix to explicitly model the spatial correlation structure across tissue locations. 

<img align="top" src="https://raw.githubusercontent.com/shangll123/workflowr_Test/main/docs/assets/main_figure.jpeg" alt="drawing" width="600"/>


## Install the Package
You can install the current version of SpatialPCA from GitHub with:
```r
library(devtools)
install_github("shangll123/SpatialPCA")
```
Please make sure you have installed the folowing R packages: 

For matrix multiplication: Matrix, RSpectra;

For spatial gene selection: SPARK;

For expression data normalization: Seurat;

For fast building large kernel matrix: parallel, pdist, tidyr, dplyr;

For result visualization: ggplot2.

## Package Tutorial
[SpatialPCA tutorial website.](http://lulushang.org/SpatialPCA_Tutorial/)

The tutorial includes main example codes for multiple spatial transcriptomics datasets (e.g. DLPFC, Slide-Seq cerebellum, Slide-Seq V2 hippocampus, Human breast tumor, and Vizgen MERFISH.)

Other analysis codes for this project can be found [here](http://lulushang.org/docs/Projects/SpatialPCA).

## Operating systems (version 1.2.0 SpatialPCA) tested on:
macOS Catalina 10.15.7

Ubuntu 18.04.5 LTS (Bionic Beaver)

## License

SpatialPCA is licensed under the GNU General Public License v3.0.

## Citation
Lulu Shang, and Xiang Zhou (2022). Spatially aware dimension reduction for spatial transcriptomics. [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.01.19.476966v1).

doi: https://doi.org/10.1101/2022.01.19.476966
