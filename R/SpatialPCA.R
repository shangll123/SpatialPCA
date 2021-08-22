#####################################################################
# Package: SpatialPCA
# Version: 1.0.1
# Date : 2021-07-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics 
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu 
#          University of Michigan, Department of Biostatistics
######################################################################
#' Each SpatialPCA object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot counts: The raw expression count matrix.
#' @slot normalized_expr: Normalized (default is SCTransform normalization from Seurat R package) expression matrix.
#' @slot project: Name of the project (for record keeping).
#' @slot covariate: The covariates in experiments (if any covariate included).
#' @slot location: Cell/spot corrdinates to compute the kernel matrix.
#' @slot kernelmat: The kernel matrix for spatial relationship between locations.
#' @slot kerneltype: The type of kernel to be used, either "gaussian", or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @slot bandwidthtype: The type of bandwidth to be used in Gaussian kernel, "SJ" for Sheather & Jones (1991) method (usually used in small size datasets), "Silverman" for Silverman's ‘rule of thumb’ method (1986)(usually used in large size datasets).
#' @slot bandwidth: The bandwidth in Gaussian kernel.
#' @slot sparseKernel: To choose if the user wants to use a sparse kernel matrix or not. It is recommended to choose sparseKernel="TRUE" when sample size is very large. 
#' @slot sparseKernel_tol: When sparseKernel=TRUE, the cut-off value when building sparse kernel matrix, any element in the kernel matrix greater than sparseKernel_tol will be kept, otherwise will be set to 0 to save memory.
#' @slot sparseKernel_ncore: When sparseKernel=TRUE, the number of CPU cores to use when building sparse kernel matrix.
#' @slot fast: Select "TRUE" to accrelerate the algorithm by performing low-rank approximation on the kernel matrix, otherwise "FALSE".
#' @slot eigenvecnum: User could optionally specify the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix.
#' @slot tau: The variance parameter in covariance matrix \Sigma, to be inferred through the algorithm.
#' @slot sigma2_0: The residual error variance, to be inferred through the algorithm.
#' @slot SpatialPCnum: The number of Spatial PCs, specified by the user, default is 20.
#' @slot tau: The variance parameter in covariance matrix for the distribution of the spatial PCs.
#' @slot W: The factor loading matrix.
#' @slot SpatialPCs: The estimated spatial PCs.
#' @slot highPCs: The estimated high resolution spatial PCs.
#' @slot highPos: The normalized positions of estimated high resolution spatial PCs.
#' @slot params: List of model parameters.
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
  params = "ANY"
) )


#' Create the SpatialPCA object with filtering and normalization step.
#' @param counts: Gene expression count matrix (matrix), m x n, m is the number of genes and n is the number of locations.
#' @param location: Spatial location matrix (matrix), n x d, d is dimensin of spatial coordinates, e.g. d=2 for locations on 2D space. 
#' @param covariate: The covariates in experiments (data.frame, if any covariate included), n x q, n is the number of locations, q is the number of covariates.
#' @param project: Project names.
#' @param min.loctions: The features (genes) detected in at least this many loctions.
#' @param min.features: The locations where at least this many features (genes) are detected.
#' @param gene.type: The type of genes to be used: "hvg" for highly variable genes; "spatial" for spatially expressed genes; "custom" for user specified genes.
#' @param gene.number: The number of top highly variable genes if gene.selection=="hvg" (all HVG genes if number is not specified); 
#' number of top spatially expressed genes if gene.selection=="spatial" (all significant spatially expressed genes if number is not specified).
#' @param customGenelist: A list of user specified genes if  gene.type=="custom".
#' @param sparkversion: In spatial gene selection, specify "spark" for small sample size data, "sparkx" for large sample size data.
#' @param numCores_spark: In spatial gene selection, specify the number of cores in SPARK package when selecting spatial genes.
#' @return Returns SpatialPCA object with filtered gene expression matrix.
#' 
#' @export
CreateSpatialPCAObject <- function(counts, location, covariate=NULL,project = "SpatialPCA", gene.type="hvg", sparkversion="spark",numCores_spark=1, gene.number=NULL,customGenelist=NULL,min.loctions = 20,  min.features=20){
  
	## check dimension
	if(ncol(counts)!=nrow(location)){
		stop("The number of cells in counts and location should be consistent! (counts -- m genes x n locations; location -- n locations x d dimension)")
	}# end fi
  
	## check data order should consistent
	if(!identical(colnames(counts), rownames(location))){
		stop("The column names of counts and row names of location should be should be matched each other! (counts -- m genes x n locations; location -- n locations x d dimension)")
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
				stop("The row names of covariate and row names of location should be should be matched each other! (covariate -- n locations x q covariates; location -- n locations x d dimension)")
			}# end fi
		
		# remove the intercept if added by user, later intercept will add automatically
			if(length(unique(covariate[,1])) == 1){
				covariate = covariate[, -1]
			}# end fi
		
		object@covariate = as.matrix(covariate)

	}# end fi

 	suppressMessages(require(Seurat))


  Seu <- CreateSeuratObject(counts = counts, project = project, min.cells = min.loctions, min.features = min.features)
	
	object@counts <- counts # store count matrix in sparse matrix
	object@location <- location
	object@project <- project


	if(!is.null(customGenelist)){
		cat(paste("## Use SCTransform function in Seurat to normalize data. \n"))
		Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
 		
 		cat(paste("## Custom gene list contains ",length(customGenelist)," genes. \n"))
 		customGenelist = as.character(customGenelist)
 		ind_match = na.omit(match(customGenelist, rownames(object@counts)))
 		cat(paste("## ",length(ind_match)," custom genes are matched with genes in the input count matrix. \n"))
 		object@normalized_expr = Seu@assays$SCT@scale.data[ind_match,]

	}else{

  if(gene.type=="hvg"){

  		cat(paste("## Use SCTransform function in Seurat to normalize data. \n"))

  	if(is.null(gene.number)){
  		
  		Seu = SCTransform(Seu, return.only.var.genes = TRUE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
  		
  		cat(paste("## Gene number is not specified, using all highly variable genes. \n"))

  		object@normalized_expr = Seu@assays$SCT@scale.data

  	}else{
  		
  		Seu = SCTransform(Seu, return.only.var.genes = TRUE, variable.features.n = gene.number,  variable.features.rv.th = 1.3)
  		
  		cat(paste("## Using top ",gene.number," highly variable genes. \n"))
  		
  		object@normalized_expr = Seu@assays$SCT@scale.data

  	}
  }else if(gene.type=="spatial"){

  	# normalize data
  	cat(paste("## Use SCTransform function in Seurat to normalize data. \n"))
  	Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
 	

  	# select spatial genes
 	if(sparkversion=="spark"){
		
cat(paste("## Use spark.test function in SPARK package to select spatially variable genes. \n"))
                suppressMessages(require(SPARK))
                count_test_spark = object@counts[match(rownames(Seu@assays$SCT@scale.data),
                  rownames(object@counts)), match(colnames(Seu@assays$SCT@scale.data),
                  colnames(object@counts))]
                location_test_spark = as.data.frame(object@location[match(colnames(Seu@assays$SCT@scale.data),
                  rownames(object@location)), ])
                spark_result <- spark(count_test_spark, location_test_spark,
                  numCores = numCores_spark)
                significant_gene_number = sum(spark_result@res_mtest$adjusted_pvalue <=
                  0.05)
                SVGnames = rownames(spark_result@res_mtest[order(spark_result@res_mtest$adjusted_pvalue),
                  ])[1:significant_gene_number]
                cat(paste("## Identified ", length(SVGnames),
                  " spatial genes through spark.test function. \n"))

 		}else if(sparkversion=="sparkx"){

 		cat(paste("## Use sparkx function in SPARK to select spatially variable genes. \n"))
 		suppressMessages(require(SPARK))
 		count_test_spark = object@counts[match(rownames(Seu@assays$SCT@scale.data), rownames(object@counts)),match(colnames(Seu@assays$SCT@scale.data), colnames(object@counts))]
 		location_test_spark = as.data.frame(object@location[match(colnames(Seu@assays$SCT@scale.data), rownames(object@location)),])
 		location_test_spark = as.matrix(location_test_spark)
 		sparkX <- sparkx(	count_test_spark, location_test_spark,numCores=numCores_spark)
 		significant_gene_number = sum(sparkX$res_mtest$adjustedPval<=0.05)
 		SVGnames = rownames(sparkX$res_mtest[order(sparkX$res_mtest$adjustedPval),])[1:significant_gene_number]
 		cat(paste("## Identified ",length(SVGnames)," spatial genes through SPARK-X function. \n"))

 	  }


 		# subset normalized data with spatial genes
  	if(is.null(gene.number)){
  		
  		cat(paste("## Gene number is not specified, using all spatially variable genes. \n"))
			object@normalized_expr = Seu@assays$SCT@scale.data[na.omit(match(SVGnames, rownames(Seu@assays$SCT@scale.data))),]

		}else {

			if(length(SVGnames) < gene.number){
				cat("The  number of significant spatial genes is less than the specified number of spatial genes! \n")
				cat(paste("## Using ",length(SVGnames)," significant spatially variable genes. \n"))
				object@normalized_expr = Seu@assays$SCT@scale.data[na.omit(match(SVGnames, rownames(Seu@assays$SCT@scale.data))),]
			}# end fi
		
			cat(paste("## Using ",gene.number," significant spatially variable genes. \n"))
			object@normalized_expr = Seu@assays$SCT@scale.data[na.omit(match(SVGnames[1:gene.number], rownames(Seu@assays$SCT@scale.data))),]
		}
  }

}

	#  location for normalized expression matrix
	object@location = object@location[match(colnames(object@normalized_expr), rownames(object@location)),]

 	  ## store as sparse matrix
	if(class(object@counts) != "dgCMatrix" ){
		object@counts <- as(object@counts, "dgCMatrix")
	}# end fi

	## covariates, i.e., confounding or batch effects
	object@params = list()
	rm(counts)
	rm(location)
	rm(Seu)


	return(object)
}# end function


#' 
spark = function(rawcount, location, numCores){
	library(SPARK)
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
                    verbose = T)
	# ress = spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
	return(spark)  
}











