########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' Calculating kernel matrix from spatial locations.
#'
#' @param object SpatialPCA object.
#' @param kerneltype The type of kernel to be used, either "gaussian", or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param bandwidthtype The type of bandwidth to be used in Gaussian kernel, "SJ" for Sheather & Jones (1991) method (usually used in small size datasets), "Silverman" for Silverman's ‘rule of thumb’ method (1986)(usually used in large size datasets).
#' @param bandwidth.set.by.user User could select their own bandwidth (a numeric value) if the recommended bandwidth doesn't work in their dataset.
#' @param sparseKernel Select "TURE" if the user wants to use a sparse kernel matrix or "FALSE" if not. It is recommended to choose sparseKernel="TRUE" when sample size is large.
#' @param sparseKernel_tol When sparseKernel=TRUE, the cut-off value when building sparse kernel matrix, any element in the kernel matrix greater than sparseKernel_tol will be kept, otherwise will be set to 0 to save memory.
#' @param sparseKernel_ncore When sparseKernel=TRUE, the number of CPU cores to build sparse kernel matrix.
#' @export
SpatialPCA_buildKernel = function(object, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL,sparseKernel=FALSE,sparseKernel_tol=1e-20,sparseKernel_ncore=1) {

	## extract the data from the slot of object, createSpatialPCAobject() function goes first
	if(length(object@counts) == 0) {
        stop("object@counts has not been set. Run CreateSpatialPCAObject() first and then retry.")
    }# end fi

	object@kerneltype = kerneltype
	object@bandwidthtype = bandwidthtype

	cat(paste("## Selected kernel type is: ", object@kerneltype," \n"))

	#************************************************************#
	#    Calculate the bandwidth for kernel matrix               #
	#************************************************************#

	#cat(paste("## Scale the expression of each gene. \n"))
	expr=object@normalized_expr
 	for(i in 1:dim(object@normalized_expr)[1]){
  		expr[i,] = scale(object@normalized_expr[i,])
	}

	if(is.null(bandwidth.set.by.user)){
		object@bandwidth = bandwidth_select(expr, method=object@bandwidthtype)
	}else{
		object@bandwidth = bandwidth.set.by.user
	}
	object@params$expr=expr

	rm(expr)

	cat(paste("## The bandwidth is: ", object@bandwidth," \n"))

  	#************************************************************#
	#    Calculate the kernel matrix with above bandwidth        #
	#************************************************************#

	location_normalized = scale(object@location)

	if(sparseKernel==FALSE){
		cat(paste("## Calculating kernel matrix\n"))
  		object@kernelmat = kernel_build(kerneltype=object@kerneltype,location=location_normalized, bandwidth=object@bandwidth)
  		object@sparseKernel=sparseKernel
  	}else if(sparseKernel==TRUE){
  		cat(paste("## Calculating sparse kernel matrix\n"))
  		object@sparseKernel=sparseKernel
  		object@sparseKernel_tol = sparseKernel_tol
		object@sparseKernel_ncore = sparseKernel_ncore
		object@kernelmat = kernel_build_sparse(kerneltype=object@kerneltype,location=location_normalized, bandwidth=object@bandwidth,tol = object@sparseKernel_tol, ncores=object@sparseKernel_ncore)

  	}

  	cat(paste("## Finished calculating kernel matrix.\n"))

	# return results
	return(object)
}# end function



#' @title Select bandwidth in Gaussian kernel.
#' @description This function selects bandwidth in Gaussian kernel.
#' @param expr A m gene by n location matrix of normalized gene expression matrix.
#' @param method The method used in bandwidth selection, "SJ" usually for small sample size data, "Silverman" usually for large sample size data.
#' @return A numeric value of calculated bandwidth.
#' @export
bandwidth_select=function (expr, method)
{
    N = dim(expr)[2]
    if (method == "SJ") {

              bw_SJ = c()
        for (i in 1:dim(expr)[1]) {
            tryCatch({
              #print(i)
            bw_SJ[i] = bw.SJ(expr[i, ], method = "dpi")
             }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
        }

        beta = median(na.omit(bw_SJ))
    }
    else if (method == "Silverman") {
        bw_Silverman = c()
        for (i in 1:dim(expr)[1]) {
            tryCatch({
            bw_Silverman[i] = bw.nrd0(expr[i, ])
            }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
        }
        beta = median(na.omit(bw_Silverman))
    }
}



#' @title Build kernel matrix.
#' @description This function calculates kernel matrix from spatial locations.
#' @param kerneltype The type of kernel to be used, either "gaussian", or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @return The kernel matrix for spatial relationship between locations.
#' @export
kernel_build = function (kerneltype = "gaussian", location, bandwidth)
{
    if (kerneltype == "gaussian") {
        K = exp(-1*as.matrix(dist(location)^2)/bandwidth)
    }
    else if (kerneltype == "cauchy") {
        K = 1/(1 + 1*as.matrix(dist(location)^2)/as.numeric(bandwidth))
    }
    else if (kerneltype == "quadratic") {
    	ED2=1*as.matrix(dist(location)^2)
        K = 1 - ED2/(ED2 + as.numeric(bandwidth))
    }
    return(K)
}




#' @title Build sparse kernel matrix.
#' @description This function calculates kernel matrix.
#' @param kerneltype The type of kernel to be used, either "gaussian", or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @param tol A numeric value of cut-off value when building sparse kernel matrix.
#' @param ncores A integer value of number of CPU cores to use when building sparse kernel matrix.
#' @return The sparse kernel matrix for spatial relationship between locations.
#'
#' @import parallel
#' @import MASS
#' @import pdist
#' @import tidyr
#'
#' @export
kernel_build_sparse = function(kerneltype,location, bandwidth,tol, ncores)
{

	# suppressMessages(require(tidyr))
	# suppressMessages(require(parallel))
	# suppressMessages(require(MASS))
	# suppressMessages(require(pdist))
	# suppressMessages(require(Matrix))

	if (kerneltype == "gaussian") {
		fx_gaussian <- function(i){
			line_i = rep(0,dim(location)[1])
			line_i[i] = 1
			line_i[-i] = exp(-(pdist(location[i,],location[-i,])@dist^2)/bandwidth)
			ind_i=which(line_i>=tol)
			return(list("ind_i"=ind_i,"ind_j"=rep(i,length(ind_i)),"val_i"=line_i[ind_i] ))
		}

		results = mclapply(1:dim(location)[1], fx_gaussian, mc.cores = ncores)
  		tib = tibble(results)  %>%  unnest_wider(results)
		K_sparse = Matrix::sparseMatrix(i =unlist(tib[[1]]), j= unlist(tib[[2]]), x= unlist(tib[[3]]),  dims = c(dim(location)[1],dim(location)[1] ))
        #K = exp(-1*as.matrix(dist(location)^2)/bandwidth)
    }
    else if (kerneltype == "cauchy") {
    	fx_cauchy <- function(i){
			line_i = rep(0,dim(location)[1])
			line_i[i] = 1
			line_i[-i] = 1/(1 + (pdist(location[i,],location[-i,])@dist^2)/as.numeric(bandwidth))
			ind_i=which(line_i>=tol)
			return(list("ind_i"=ind_i,"ind_j"=rep(i,length(ind_i)),"val_i"=line_i[ind_i] ))
		}

		results = mclapply(1:dim(location)[1], fx_cauchy, mc.cores = ncores)
  		tib = tibble(results)  %>%  unnest_wider(results)
		K_sparse = Matrix::sparseMatrix(i =unlist(tib[[1]]), j= unlist(tib[[2]]), x= unlist(tib[[3]]),  dims = c(dim(location)[1],dim(location)[1] ))
        # K = 1/(1 + 1*as.matrix(dist(location)^2)/as.numeric(bandwidth))
    }
    else if (kerneltype == "quadratic") {

    	fx_quadratic <- function(i){
			line_i = rep(0,dim(location)[1])
			line_i[i] = 1
			ED2=pdist(location[i,],location[-i,])@dist^2
			line_i[-i] = 1 - ED2/(ED2 + as.numeric(bandwidth))
			ind_i=which(line_i>=tol)
			return(list("ind_i"=ind_i,"ind_j"=rep(i,length(ind_i)),"val_i"=line_i[ind_i] ))
		}

		results = mclapply(1:dim(location)[1], fx_quadratic, mc.cores = ncores)
  		tib = tibble(results)  %>%  unnest_wider(results)
		K_sparse = sparseMatrix(i =unlist(tib[[1]]), j= unlist(tib[[2]]), x= unlist(tib[[3]]),  dims = c(dim(location)[1],dim(location)[1] ))
       # ED2=1*as.matrix(dist(location)^2)
       # K = 1 - ED2/(ED2 + as.numeric(bandwidth))

    }
    return(K_sparse)

}



















