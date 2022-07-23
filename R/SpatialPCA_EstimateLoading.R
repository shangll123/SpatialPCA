########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' Calculate loading matrix.
#'
#' @param object SpatialPCA object.
#' @param maxiter Maximum iteration number. Default is 300.
#' @param initial_tau Initial value of tau. Default is 1. Because we need tau to be positive, we calculate exp(log(tau)) during iterations.
#' @param fast Select "TRUE" if the user wants to use low-rank approximation on the kernel matrix to accelerate the algorithm, otherwise select "FALSE".
#' @param eigenvecnum When fast=TRUE, eigenvecnum is the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix.
#' The default is NULL, if specified, it is recommended to use eigenvecnum=20 when sample size is large (e.g. >5,000). When sample size is small, eigenvecnum is suggested to explain at least 90% variance.
#' @param SpatialPCnum Number of spatial PCs.
#' @return Returns SpatialPCA object with estimated loading matrix W.
#'
#' @import RSpectra
#'
#' @export
#'
SpatialPCA_EstimateLoading = function(object, maxiter=300,initial_tau=1,fast=FALSE,eigenvecnum=NULL,SpatialPCnum=20){

      suppressMessages(require(RSpectra))
      set.seed(1234)
      param_ini=log(initial_tau)
      object@SpatialPCnum = SpatialPCnum
      object@fast = fast
      object@params$X = scale(object@location)
      object@params$n = dim(object@params$X)[1]
      object@params$p=dim(object@params$X)[2]

  if(is.null(object@covariate)){
      object@params$H = matrix(1, dim(object@params$X)[1],1)
      HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
      HH = object@params$H%*%HH_inv%*%t(object@params$H)
      object@params$M=diag(object@params$n)-HH
      # Y=expr
      object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
      object@params$YM = object@params$expr%*%object@params$M
      object@params$q=1
  }else{
      object@params$q = dim(object@covariate)[2]+1
      object@params$H = matrix(0, object@params$n,object@params$q)
      object@params$H[,1]=1
      object@params$H[,2:object@params$q] = object@covariate
      HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
      HH=object@params$H%*%HH_inv%*%t(object@params$H)
      object@params$M=diag(object@params$n)-HH
      #Y=expr
      object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
      object@params$YM = object@params$expr%*%object@params$M
  }


  if(fast==FALSE){
      object@fast=fast
      print("Eigen decomposition on kernel matrix!")
      eigen_res = eigen(object@kernelmat)
      object@params$delta = eigen_res$values
      object@params$U = eigen_res$vectors
      print("Using all eigenvectors and eigenvalues in the Kernel matrix!")
  }else{
      object@fast=fast
      if(!is.null(eigenvecnum)){
        print("Eigen decomposition on kernel matrix!")
        object@eigenvecnum=eigenvecnum
        if(object@sparseKernel==TRUE){
          eigen_res = eigs_sym(object@kernelmat, k=object@eigenvecnum)
          object@params$delta = eigen_res$values
          object@params$U = eigen_res$vectors
          }else{
          eigen_res = eigs_sym(object@kernelmat, k=object@eigenvecnum, which = "LM")
          object@params$delta = eigen_res$values
          object@params$U = eigen_res$vectors
        }

        print("Low rank approximation!")
        print(paste0("Using user selected top ",object@eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
      }else if(object@params$n>5000){
        print("Eigen decomposition on kernel matrix!")
        if(object@sparseKernel==TRUE){
          eigen_res = eigs_sym(object@kernelmat, k=20)
          object@params$delta = eigen_res$values
          object@params$U = eigen_res$vectors
          }else{
          eigen_res = eigs_sym(object@kernelmat, k=20, which = "LM")
          object@params$delta = eigen_res$values
          object@params$U = eigen_res$vectors
        }
        print("Low rank approximation!")
        print("Large sample, using top 20 eigenvectors and eigenvalues in the Kernel matrix!")
        }else{
        eigen_res = eigen(object@kernelmat)
        delta_all = eigen_res$values
        U_all = eigen_res$vectors
        ind = which(cumsum(delta_all/length(delta_all))>0.9)[1]
        print("Low rank approximation!")
        print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
        object@params$delta = delta_all[1:ind]
        object@params$U = U_all[,1:ind]
        rm(U_all)
        }
    }


    object@params$MYt = object@params$M %*% t(object@params$expr)
    object@params$YMMYt = object@params$YM %*% object@params$MYt
    object@params$YMU = object@params$YM %*% object@params$U
    object@params$Xt = t(object@params$H)
    object@params$XtU = object@params$Xt %*% object@params$U
    object@params$Ut = t(object@params$U)
    object@params$UtX = object@params$Ut %*% object@params$H
    object@params$YMX = object@params$YM %*% object@params$H
    object@params$UtU = object@params$Ut %*% object@params$U
    object@params$XtX = object@params$Xt %*% object@params$H
    object@params$SpatialPCnum = SpatialPCnum


    optim_result =try(optim(param_ini, SpatialPCA_estimate_parameter,params=object@params,control = list(maxit = maxiter), lower = -10, upper = 10,method="Brent"),silent=T)

    object@tau = exp(optim_result$par)
    k = dim(object@params$expr)[1]
    n = dim(object@params$expr)[2]
    q=object@params$q
    tauD_UtU_inv = solve(object@tau*diag(object@params$delta) + object@params$UtU, tol = 1e-40)
    YMU_tauD_UtU_inv_Ut = object@params$YMU %*% tauD_UtU_inv %*% object@params$Ut
    YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% object@params$H
    XtU_inv_UtX = object@params$XtU %*% tauD_UtU_inv %*% object@params$UtX
    left = object@params$YMX - YMU_tauD_UtU_inv_UtX
    right = t(left)
    middle = solve(-XtU_inv_UtX, tol = 1e-40)
    G_each = object@params$YMMYt - YMU_tauD_UtU_inv_Ut %*% object@params$MYt - left %*% middle %*% right
    object@W = eigs_sym(G_each, k=SpatialPCnum, which = "LM")$vectors
    object@sigma2_0 = as.numeric((object@params$tr_YMY+F_funct_sameG(object@W,G_each))/(k*(n-q)))

    rm(eigen_res)
    rm(tauD_UtU_inv)
    rm(YMU_tauD_UtU_inv_Ut)
    rm(YMU_tauD_UtU_inv_UtX)
    rm(XtU_inv_UtX)
    rm(left)
    rm(right)
    rm(middle)
    rm(G_each)
    gc()

  return(object)
}


#' @import RSpectra
SpatialPCA_estimate_parameter = function(param_ini, params){
    # suppressMessages(require(RSpectra))
    set.seed(1234)
    tau=exp(param_ini[1])
    k = dim(params$expr)[1]
    n = dim(params$expr)[2]
    q=params$q
    PCnum=params$SpatialPCnum
    tauD_UtU_inv = solve(tau*diag(params$delta) + params$UtU, tol = 1e-40)
    YMU_tauD_UtU_inv_Ut = params$YMU %*% tauD_UtU_inv %*% params$Ut
    YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% params$H
    XtU_inv_UtX = params$XtU %*% tauD_UtU_inv %*% params$UtX
    left = params$YMX - YMU_tauD_UtU_inv_UtX
    right = t(left)
    middle = solve(-XtU_inv_UtX, tol = 1e-40)
    G_each = params$YMMYt - YMU_tauD_UtU_inv_Ut %*% params$MYt - left %*% middle %*% right
    log_det_tauK_I = determinant(1/tau*diag(1/params$delta)+ params$UtU, logarithm=TRUE)$modulus[1] + determinant(tau*diag(params$delta), logarithm=TRUE)$modulus[1]
    Xt_invmiddle_X = params$XtX - params$XtU %*% solve(params$UtU + 1/tau *diag( 1/params$delta) , tol = 1e-40) %*% params$UtX
    log_det_Xt_inv_X = determinant(Xt_invmiddle_X, logarithm=TRUE)$modulus[1]
    sum_det=0
    sum_det=sum_det+(0.5*log_det_tauK_I+0.5*log_det_Xt_inv_X  )*PCnum

    rm(tauD_UtU_inv)
    rm(YMU_tauD_UtU_inv_Ut)
    rm(YMU_tauD_UtU_inv_UtX)
    rm(XtU_inv_UtX)
    rm(left)
    rm(middle)
    rm(right)
    rm(Xt_invmiddle_X)
    gc()

    W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
    -(-sum_det -(k*(n-q))/2*log(params$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}


F_funct_sameG = function(X,G){ # G is a matrix
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G%*%X[,i]
  }
  -return_val
}




