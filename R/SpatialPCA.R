#####################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics 
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu 
#          University of Michigan, Department of Biostatistics
######################################################################

#' Each SpatialPCA object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot counts The raw expression count matrix. Rows are genes, columns are spots/cells.
#' @slot normalized_expr Normalized (by default we use SCTransform normalization in Seurat R package) expression matrix.
#' @slot project Name of the project (for record keeping).
#' @slot covariate The covariates in experiments (if any covariate included).
#' @slot location Cell/spot spatial coordinates to compute the kernel matrix.
#' @slot kernelmat The kernel matrix for spatial relationship between locations.
#' @slot kerneltype The type of kernel to be used, either "gaussian" for gaussian kernel, or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @slot bandwidthtype The type of bandwidth to be used in Gaussian kernel, "SJ" for Sheather & Jones (1991) method (usually used in small sample size datasets), "Silverman" for Silverman's ‘rule of thumb’ method (1986)(usually used in large sample size datasets).
#' @slot bandwidth The bandwidth in Gaussian kernel, users can also specify their preferred bandwidth.
#' @slot sparseKernel To choose if the user wants to use a sparse kernel matrix or not. It is recommended to choose sparseKernel="TRUE" when sample size is large and you want to speed up the calculation. 
#' @slot sparseKernel_tol When sparseKernel=TRUE, sparseKernel_tol is the cut-off value when building sparse kernel matrix, any element in the kernel matrix greater than sparseKernel_tol will be kept, otherwise will be set to 0 to save memory.
#' @slot sparseKernel_ncore When sparseKernel=TRUE, sparseKernel_ncore is the number of CPU cores to use when building the sparse kernel matrix.
#' @slot fast Select "TRUE" to accrelerate the algorithm by performing low-rank approximation on the kernel matrix, otherwise "FALSE" for calculation without low-rank approximation on the kernel matrix.
#' @slot eigenvecnum When fast=TRUE, the user can optionally specify the number of top eigenvectors and eigenvalues to be used in low-rank approximation when performing eigen decomposition on the kernel matrix.
#' @slot tau The variance parameter in covariance matrix for the spatial PCs, to be inferred through the algorithm.
#' @slot sigma2_0 The residual error variance, to be inferred through the algorithm.
#' @slot SpatialPCnum The number of Spatial PCs, specified by the user, default is 20.
#' @slot W The factor loading matrix.
#' @slot SpatialPCs The estimated spatial PCs.
#' @slot highPCs The estimated high resolution spatial PCs, if needed.
#' @slot highPos The scaled locations of estimated high resolution spatial PCs, if needed.
#' @slot expr_pred The predicted gene expression on new locations when highPCs and highPos are avaliable.
#' @slot params List of model parameters.
#' @export

setClass("SpatialPCA", slots=list(
  counts = "ANY",
  normalized_expr = "ANY",
  project = "character",
  covariate = "ANY",
  location = "matrix", 
  kernelmat = "ANY",
  kerneltype = "character",
  bandwidthtype = "character",
  bandwidth = "numeric",
  sparseKernel="logical",
  sparseKernel_tol = "numeric",
  sparseKernel_ncore = "numeric",
  fast = "logical",
  eigenvecnum = "numeric",
  SpatialPCnum = "numeric",
  tau = "numeric",
  sigma2_0 = "numeric",
  W = "ANY",
  SpatialPCs = "ANY",
  highPCs = "ANY",
  highPos = "ANY",
  expr_pred="ANY",
  params = "ANY"
) )


#' Create the SpatialPCA object with filtering and normalization step.
#' @param counts Gene expression count matrix (matrix), the dimension is m x n, where m is the number of genes and n is the number of locations.
#' @param location Spatial location matrix (matrix), the dimension is n x d, n is the number of locations, d is dimensin of spatial coordinates, e.g. d=2 for locations on 2D space. The rownames of locations and the colnames of count matrix should be matched.
#' @param covariate The covariates in experiments (matrix, if any covariate included), n x q, n is the number of locations, q is the number of covariates. The rownames of covariates and the rownames of locations should be matched.
#' @param project Name of the project (for record keeping).
#' @param min.loctions The features (genes) detected in at least min.loctions number of loctions, default is 20.
#' @param min.features The locations where at least min.features number of features (genes) are detected, default is 20.
#' @param gene.type The type of genes to be used: "spatial" for spatially expressed genes; "hvg" for highly variable genes; "custom" for user specified genes, default is "spatial".
#' @param gene.number The number of top highly variable genes if gene.selection=="hvg" (use all HVG genes if this number is not specified); 
#' number of top spatially expressed genes if gene.selection=="spatial" (use all significant spatially expressed genes if this number is not specified).
#' @param customGenelist A list of user specified genes if  gene.type=="custom".
#' @param sparkversion In spatial gene selection, specify "spark" for small sample size data for higher detection power of spatial genes, "sparkx" for large sample size data for saving time and memory.
#' @param numCores_spark If gene.type="spatial", specify the number of CPU cores in SPARK package to use when selecting spatial genes.
#' @return Returns SpatialPCA object, with filtered and normalized gene expression matrix and corresponding location matrix.
#' 
#' @import Seurat
#' @import SPARK
#' 
#' @examples 
#' 
#' 
#' 
#' @export
CreateSpatialPCAObject <- function(counts, location, covariate=NULL,project = "SpatialPCA", gene.type="spatial", sparkversion="spark",numCores_spark=1, gene.number=3000,customGenelist=NULL,min.loctions = 20,  min.features=20){
  
  #suppressMessages(require(Seurat))
	#suppressMessages(require(SPARK))

	## check dimension
	if(ncol(counts)!=nrow(location)){
		stop("The number of cells in counts and location should be consistent (counts -- m genes x n locations; location -- n locations x d dimension).")
	}# end fi
  
	## check data order should consistent
	if(!identical(colnames(counts), rownames(location))){
		stop("The column names of counts and row names of location should be should be matched (counts -- m genes x n locations; location -- n locations x d dimension).")
	}# end fi
 

	## inheriting
	object <- new(
		Class = "SpatialPCA",
		counts = counts,
		location = location,
		project = project
	)
  
    if(!is.null(covariate)){
		  	## check data order should consistent
			if(!identical(rownames(covariate), rownames(location))){
				stop("The row names of covariate and row names of location should be should be matched (covariate -- n locations x q covariates; location -- n locations x d dimension).")
			}# end fi
			
			q=dim(covariate)[2]
			n_covariate=dim(covariate)[1]
		# remove the intercept if added by user, later intercept will add automatically
			if(length(unique(covariate[,1])) == 1){
				covariate = covariate[, -1]
				q=q-1
			}# end fi
		
		object@covariate = as.matrix(covariate,n_covariate,q)

	}# end fi

 

  Seu <- CreateSeuratObject(counts = counts, project = project, min.cells = min.loctions, min.features = min.features)
	
	object@counts <- counts # store count matrix in sparse matrix
	object@location <- location
	object@project <- project

	rm(counts) # to save memory
	rm(location)

	if(!is.null(customGenelist)){ # if user specified customGenelist

								cat(paste("## Use SCTransform function in Seurat to normalize data. \n"))
								Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3,verbose = FALSE)
						 		cat(paste("## Custom gene list contains ",length(customGenelist)," genes. \n"))
						 		customGenelist = as.character(customGenelist)
						 		ind_match = na.omit(match(customGenelist, rownames(object@counts)))
						 		cat(paste("## In total ",length(ind_match)," custom genes are matched with genes in the count matrix. \n"))
						 		object@normalized_expr = Seu@assays$SCT@scale.data[ind_match,]
						 		cat(paste("## Use ",length(ind_match)," custom genes for analysis. \n"))

	}else{		# if user didn't specify customGenelist

  	if(gene.type=="hvg"){

								cat(paste("## Use SCTransform function in Seurat to normalize data. \n"))

								if(is.null(gene.number)){

								Seu = SCTransform(Seu, return.only.var.genes = TRUE, variable.features.n = NULL,  variable.features.rv.th = 1.3,verbose = FALSE)
								object@normalized_expr = Seu@assays$SCT@scale.data
								gene.number = dim(Seu@assays$SCT@scale.data)[1]
								cat(paste("## Gene number is not specified, using all ",gene.number," highly variable genes. \n"))  		

								}else{

								Seu = SCTransform(Seu, return.only.var.genes = TRUE, variable.features.n = NULL,  variable.features.rv.th = 1.3,verbose = FALSE)
								hvg_gene_num = dim(Seu@assays$SCT@scale.data)[1]

								if( gene.number < hvg_gene_num ){

								object@normalized_expr = Seu@assays$SCT@scale.data[1:gene.number,]
								
								cat(paste("## Using top ",gene.number," highly variable genes. \n"))

								}else{

								object@normalized_expr = Seu@assays$SCT@scale.data
								cat("The  number of highly variable genes is less than the specified number of genes. \n")
								cat(paste("## Using ",hvg_gene_num," highly variable genes for analysis. \n"))

								}

								}
  	}else if(gene.type=="spatial"){

						  	# normalize data
						  	cat(paste("## Use SCTransform function in Seurat to normalize data. \n"))
						  	Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3,verbose = FALSE)
						 	
				  	# select spatial genes
				 	if(sparkversion=="spark"){
												cat(paste("## Use spark.test function in SPARK package to select spatially variable genes. \n"))
				                #suppressMessages(require(SPARK))
				                count_test_spark = object@counts[na.omit(match(rownames(Seu@assays$SCT@scale.data), rownames(object@counts))), na.omit(match(colnames(Seu@assays$SCT@scale.data),colnames(object@counts)))]
				                location_test_spark = as.data.frame(object@location[match(colnames(Seu@assays$SCT@scale.data), rownames(object@location)), ])
				                spark_result <- spark(count_test_spark, location_test_spark,numCores = numCores_spark)
				                significant_gene_number = sum(spark_result@res_mtest$adjusted_pvalue <= 0.05)
				                SVGnames = rownames(spark_result@res_mtest[order(spark_result@res_mtest$adjusted_pvalue),])[1:significant_gene_number]
				                cat(paste("## Identified ", length(SVGnames)," spatial genes through spark.test function. \n"))
					}else if(sparkversion=="sparkx"){

												cat(paste("## Use sparkx function in SPARK to select spatially variable genes. \n"))
												count_test_spark = object@counts[na.omit(match(rownames(Seu@assays$SCT@scale.data), rownames(object@counts))),na.omit(match(colnames(Seu@assays$SCT@scale.data), colnames(object@counts)))]
												location_test_spark = as.data.frame(object@location[match(colnames(Seu@assays$SCT@scale.data), rownames(object@location)),])
												location_test_spark = as.matrix(location_test_spark)
												sparkX <- sparkx(	count_test_spark, location_test_spark,numCores=numCores_spark)
												significant_gene_number = sum(sparkX$res_mtest$adjustedPval<=0.05)
												SVGnames = rownames(sparkX$res_mtest[order(sparkX$res_mtest$adjustedPval),])[1:significant_gene_number]
												cat(paste("## Identified ",length(SVGnames)," spatial genes through SPARK-X function. \n"))
					}

					# subset normalized data with spatial genes
					if(is.null(gene.number)){
						
						object@normalized_expr = Seu@assays$SCT@scale.data[na.omit(match(SVGnames, rownames(Seu@assays$SCT@scale.data))),]
						cat(paste("## Gene number is not specified, we use all ",gene.number," spatially variable genes. \n")) 
					
					}else {

						if(length(SVGnames) < gene.number){
							cat("The  number of significant spatial genes is less than the specified number of spatial genes. \n")
							cat(paste("## Using ",length(SVGnames)," significant spatially variable genes. \n"))
							object@normalized_expr = Seu@assays$SCT@scale.data[na.omit(match(SVGnames, rownames(Seu@assays$SCT@scale.data))),]
						}else{
							cat(paste("## Using top ",gene.number," significant spatially variable genes. \n"))
							object@normalized_expr = Seu@assays$SCT@scale.data[na.omit(match(SVGnames[1:gene.number], rownames(Seu@assays$SCT@scale.data))),]
						}
					}

}

}

					#  location for normalized expression matrix
					object@location = object@location[match(colnames(object@normalized_expr), rownames(object@location)),]

					#  covariates, i.e., confounding or batch effects
					if(!is.null(covariate)){
						object@covariate = object@covariate[match(colnames(object@normalized_expr), rownames(object@location)),1:q]
						object@covariate = as.matrix(object@covariate,dim(object@normalized_expr)[2],q )
					}

				 	## store count matrix as a sparse matrix
					if(class(object@counts) != "dgCMatrix" ){
						object@counts <- as(object@counts, "dgCMatrix")
					}# end fi

					
					object@params = list()

					rm(Seu)

					return(object)
}# end function


#' @import SPARK
spark = function(rawcount, location, numCores){
	# library(SPARK)
	location = as.data.frame(location)
	rownames(location) = colnames(rawcount)

	spark <- CreateSPARKObject(counts=rawcount, location=location,percentage = 0.1, min_total_counts = 10)
	spark@lib_size <- apply(rawcount, 2, sum)
	spark <- spark.vc(spark, 
          covariates = NULL, 
                    lib_size = spark@lib_size, 
                    num_core = numCores,
                    verbose = F,
                    fit.model="gaussian")
	spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = F)
	return(spark)  
}











