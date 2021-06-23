# **SpatialPCA** 
Welcome to SpatialPCA, a spatially aware dimension reduction method that aims to infer a low dimensional representation of the gene expression data in spatial transcriptomics. SpatialPCA builds upon the probabilistic version of PCA, incorporates localization information as additional input, and uses a kernel matrix to explicitly model the spatial correlation structure across tissue locations. 

<img align="top" src="https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/SpatialPCA_Figure1.png" alt="drawing" width="600"/>

Detailed tutorial could be found here: [SpatialPCA tutorial](http://lulushang.org/SpatialPCA.html)

All analysis codes could be found here: [Analysis codes](http://lulushang.org/docs/Projects/SpatialPCA)

## Install the Package
You can install the current version of SpatialPCA from GitHub with:
```r
library(devtools)
install_github("shangll123/SpatialPCA")
```

## Dependencies
R version >= 3.5.0.
R packages: ggplot2, RSpectra, bluster, slingshot.

## Installation time: 
Less than one minute, after installing dependent packages. 

## Runtime: 
The example dataset provided ([SpatialPCA tutorial](http://lulushang.org/SpatialPCA.html)) can be run in less than 2 minutes on a normal desktop computer. The extraction of spatial PCs takes approximately 27 seconds on a osmFISH dataset (33 genes and 5,275 locations) on the server using a single thread on an Intel(R) Xeon(R) Gold 6138 CPU @ 2.00GHz processor. The extraction of spatial PCs takes approximately 68min on a slide-seq dataset (787 genes and 20,982 locations) on the server using a single thread on the same processor.

## Operating systems (version 1.0.0 SpatialPCA) tested on:
macOS Catalina 10.15.7
Ubuntu 18.04.5 LTS (Bionic Beaver)

## License

SpatialPCA is licensed under the GNU General Public License v3.0.
