########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics 
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu 
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' Calculating Spatial PCs (latent factor matrix Z).
#' @param object SpatialPCA object.
#' @param fast Select fast=TRUE if the user wants to use low-rank approximation on the kernel matrix to calculate the spatial PCs, otherwise select FALSE. 
#' @param eigenvecnum When fast=TRUE, eigenvecnum is the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix. 
#' The default is NULL, if specified, it is recommended that these top eigen values explain >=90% of the variance. 
#' In estimating spatial PCs, we need larger number of eigenvectors in kernel matrix for more accurate estimation.
#' @return Returns SpatialPCA object with estimated Spatial PCs.
#' 
#' @import RSpectra
#' 
#' @export
SpatialPCA_SpatialPCs= function(object,fast=FALSE,eigenvecnum=NULL){

	# suppressMessages(require(RSpectra))

    n = object@params$n
    PCnum = object@SpatialPCnum
    Z_hat = matrix(0, PCnum, n)
    tau = object@tau
    W_hat = object@W

if(fast==FALSE){
        U=object@params$U
        delta=object@params$delta
}else if(fast==TRUE){
    
    if(!is.null(eigenvecnum)){
        print(paste0("Low rank approximation!"))
        print(paste0("Using user selected top ",eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
        EIGEN = eigs_sym(object@kernelmat, k=eigenvecnum, which = "LM")
        U=EIGEN$vectors
        delta=EIGEN$values

    }else if(n>5000){
        fast_eigen_num = ceiling(n*0.1)
        print(paste0("Low rank approximation!"))
        print("Large sample, using top 10% sample size of eigenvectors and eigenvalues in the Kernel matrix!")
        EIGEN = eigs_sym(object@kernelmat, k=fast_eigen_num, which = "LM")
        U=EIGEN$vectors
        delta=EIGEN$values
    }else{
        U=object@params$U
        delta=object@params$delta
        ind=length(delta)
        print(paste0("Low rank approximation!"))
        print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
    }
}
    object@params$U=U
    object@params$delta=delta

    W_hat_t = t(W_hat)
    WtYM = W_hat_t%*% object@params$YM
    WtYMK = WtYM %*% object@kernelmat
    WtYMU = WtYM %*% object@params$U
    Ut=t(object@params$U)
    UtM = Ut %*% object@params$M
    UtMK = UtM %*% object@kernelmat
    UtMU = UtM %*% object@params$U
    middle_inv = solve(1/tau * diag(1/delta) + UtMU, tol = 1e-40)
    
    object@SpatialPCs = tau*WtYMK - tau*WtYMU %*% middle_inv %*% UtMK

    rm(W_hat_t)
    rm(WtYM)
    rm(WtYMK)
    rm(WtYMU)
    rm(Ut)
    rm(UtM)
    rm(UtMK)
    rm(UtMU)
    rm(middle_inv)
    gc()


    return(object)
}








