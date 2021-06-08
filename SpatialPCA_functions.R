

#' @title Scale expression values of each gene
#' @description This function scales the expression of each input gene. 
#' @param expr: A m gene by n location of normalized gene expression matrix.
#' @return A m gene by n location of gene expression matrix with each row of gene scaled
#' @export
scale_expr = function(expr){
print(paste0("Input expression data: ",dim(expr)[1]," genes on ",dim(expr)[2]," locations."))
  for(i in 1:dim(expr)[1]){
  expr[i,] = scale(expr[i,])
}
  expr
}


#' @title Calculate and prepare data.
#' @description This function calculates needed information for SpatialPCA.
#' @param expr: A m gene by n location of gene expression matrix.
#' @param info: A n cell by k dimension of location matrix. n is cell number, k is dimension (k=2 when location is 2D; k=3 when location is 3D)
#' @param covariate: A covariate by n location matrix. If there is no covariate, the default is "NA". 
#' @param kerneltype: An character string. It is the type of kernel to be used in SpatialPCA. Default is "gaussian", other options include "cauchy" for cauchy kernel and "quadratic" for rational quadratic kernel.
#' @param bandwidthtype: A character string representing method used in bandwidth selection: "SJ" for Sheather & Jones (1991) method, 
#' "Silverman" for Silverman's ‘rule of thumb’ method (1986), "Scott" for Scott (1992) method.
#' @param fast: A logical value, select "TRUE" to accrelerate the algorithm by performing low-rank approximation on the kernel matrix, otherwise "FALSE". 
#' @param eigenvecnum: An integer, the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix. 
#' @return A list of data matrices.
#' \item{K}{Kernel matrix.}
#' \item{ED}{ED distance matrix.}
#' \item{ED2}{ED2 distance square matrix.}
#' \item{n}{Sample size.}
#' \item{Y}{Expression matrix.}
#' \item{q}{Number of covariates, if no input covariate is avaliable then q=1 because we added an intercept in the function.}
#' \item{bandwidth}{Bandwidth in the kernel matrix.}
#' \item{delta}{Eigenvalues of the kernel matrix.}
#' \item{U}{Eigenvectors of the kernel matrix.}
#' \item{kerneltype}{The type of kernel function used.}
#' \item{bandwidthtype}{The type of bandwidth selection method used.}
#' \item{location}{A n cell by k dimension of location matrix.}
#' \item{M}{M matrix.}
#' \item{tr_YMY}{tr_YMY, a scalar.}
#' \item{YM}{YM matrix.}
#' \item{MYt}{MYt matrix.}
#' \item{YMMYt}{YMMYt matrix.}
#' \item{YMU}{YMU matrix.}
#' \item{H}{X matrix.}
#' \item{Xt}{Xt matrix.}
#' \item{XtU}{XtU matrix.}
#' \item{Ut}{Ut matrix.}
#' \item{UtX}{UtX matrix.}
#' \item{YMX}{YMX matrix.}
#' \item{UtU}{UtU matrix.}
#' \item{XtX}{XtX matrix.}
#' @export
data_prepare_func = function(expr,info,covariate=NA,  kerneltype = "gaussian",bandwidthtype="SJ",fast = FALSE,eigenvecnum=NA){

  if (ncol(expr) != nrow(info)) {
    stop("ERROR - expression samples size doesn't agree with the number of locations")
  }

suppressMessages(require(RSpectra))

expr=scale_expr(expr)
X = scale(info)
n = dim(X)[1]
p=dim(X)[2] 
if(sum(is.na(covariate))==1){
  H = matrix(1, dim(X)[1],1)
  HH_inv=solve(t(H)%*%H) 
  HH=H%*%HH_inv%*%t(H)
  M=diag(n)-HH
  Y=expr
  tr_YMY=sum(diag(Y%*%M%*%t(Y)))
  ED = as.matrix(dist(scale(info[ ,1:2])))
  ED2 = ED^2
  YM = Y%*%M

  q=1
}else{
  if (ncol(covariate) != nrow(info)) {
    stop("ERROR - covariate samples size doesn't agree with the number of locations")
  }

  q = dim(covariate)[1]+1
  H = matrix(0, dim(X)[1],q)
  H[,1]=1
  H[,2:q] = covariate
  HH_inv=solve(t(H)%*%H) 
  HH=H%*%HH_inv%*%t(H)
  M=diag(n)-HH
  Y=expr
  tr_YMY=sum(diag(Y%*%M%*%t(Y)))
  ED = as.matrix(dist(scale(info[ ,1:2])))
  ED2 = ED^2
  YM = Y%*%M

}


print("Kernel matrix!")
  bandwidth = bandwidth_select(expr, info,method=bandwidthtype)
  K=kernel_build(kernelpara=kerneltype, ED2=ED2,bandwidth) 
  
if(fast==FALSE){
  print("Eigen decomposition on kernel matrix!")
  eigen_res = eigen(K)
  delta = eigen_res$values
  U = eigen_res$vectors
  print("Using all eigenvectors and eigenvalues in the Kernel matrix!")
}else{
    if(is.na(eigenvecnum)==FALSE){
        print("Eigen decomposition on kernel matrix!")
        eigen_res = eigs_sym(K, k=eigenvecnum, which = "LM")
        delta = eigen_res$values
        U = eigen_res$vectors
        print("Low rank approximation!")
        print(paste0("Using user selected top ",eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
    }else if(n>5000){
        print("Eigen decomposition on kernel matrix!")
        eigen_res = eigs_sym(K, k=20, which = "LM")
        delta = eigen_res$values
        U = eigen_res$vectors
        print("Low rank approximation!")
        print("Large sample, using top 20 eigenvectors and eigenvalues in the Kernel matrix!")
      }else{
        eigen_res = eigen(K)
        delta_all = eigen_res$values
        U_all = eigen_res$vectors
        ind = which(cumsum(delta_all/length(delta_all))>0.9)[1]
        print("Low rank approximation!")
        print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
        delta = delta_all[1:ind]
        U = U_all[,1:ind]
      }
    }

# data prepare
print("Preparing data for later calculation!")
    YM = Y %*% M
    MYt = M%*%t(Y)
    YMMYt = YM %*% MYt
    YMU = YM %*% U
    Xt = t(H)
    XtU = Xt %*% U
    Ut = t(U)
    UtX = Ut %*% H
    YMX = YM %*% H
    UtU = Ut %*% U
    XtX = Xt %*% H

  res<-list("K"=K,"ED"=ED, "ED2"=ED2,  "n"=n,"Y"=expr,
    "q"=q,"bandwidth"=bandwidth,"delta"=delta,"U"=U,kerneltype = "gaussian",
    bandwidthtype=bandwidthtype,location=info, "M"=M, "tr_YMY"=tr_YMY,"YM"=YM, 
    "MYt"=MYt,"YMMYt"=YMMYt,"YMU"=YMU,"H"=H,"Xt" = Xt, 
    "XtU" = XtU, "Ut" = Ut, "UtX" = UtX, "YMX"=YMX, "UtU" = UtU,
    "XtX" = XtX)
  return(res)

}



#' @title Select bandwidth in Gaussian kernel.
#' @description This function selects bandwidth in Gaussian kernel.
#' @param expr: A m gene by n location of normalized and scaled gene expression matrix.
#' @param info: A n cell by k dimension of location matrix. 
#' @param method: Method used in bandwidth selection, "SJ" for Sheather & Jones (1991) method, 
#' "Silverman" for Silverman's ‘rule of thumb’ method (1986), "Scott" for Scott (1992) method.
#' @return A numeric value of calculated bandwidth.
#' @export

bandwidth_select=function (expr, info, method) 
{
    N = dim(info)[1]
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
    else if (method == "Scott") {
        bw_Scott = c()
        for (i in 1:dim(expr)[1]) {
              tryCatch({ 
            bw_Scott[i] = bw.nrd(expr[i, ])
             }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
        }
        beta = median(na.omit(bw_Scott))
    }
    return(beta)
}




#' @title Build kernels.
#' @description This function builds kernel matrix.
#' @param kernelpara: An character string that represents different kernels. Default is "gaussian", other options include "cauchy" for cauchy kernel and "quadratic" for rational quadratic kernel.
#' @param ED2: An n by n distance squared matrix.
#' @param beta: A numeric value of bandwidth.
#' @return An n by n kernel matrix.
#' @export
kernel_build = function(kernelpara="gaussian", ED2,beta){
  if(kernelpara == "gaussian"){
    K = exp(-1.0*ED2/as.numeric(beta))
}else if(kernelpara == "cauchy"){
    K = 1.0/(1.0+ED2/as.numeric(beta))
}else if(kernelpara == "quadratic"){
    K=  1.0-ED2/(ED2+as.numeric(beta))
}
  return(K)
}



#' @title Estimation of parameters.
#' @description This function calculate -log likelihood value.
#' @param param_ini: Input log(tau) value. Because we need tau to be positive, we calculate exp(log(tau)) inside of this function.
#' @param dat_input: Object created in the data_prepare_func function.
#' @param PCnum: Number of spatial PCs.
#' @return A numeric value of -log likelihood.
#' @export
SpatialPCA_estimate = function(param_ini,dat_input,PCnum=PCnum){
   suppressMessages(require(RSpectra))
    set.seed(1234)
    param = param_ini
    #print(param)
    tau=exp(param[1])
    k = dim(dat_input$Y)[1]
    n = dim(dat_input$Y)[2]
    q=dat_input$q

    tauD_UtU_inv = solve(tau*diag(dat_input$delta) + dat_input$UtU)
    YMU_tauD_UtU_inv_Ut = dat_input$YMU %*% tauD_UtU_inv %*% dat_input$Ut
    YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% dat_input$H
    XtU_inv_UtX = dat_input$XtU %*% tauD_UtU_inv %*% dat_input$UtX
    left = dat_input$YMX - YMU_tauD_UtU_inv_UtX
    right = t(left)
    middle = solve(-XtU_inv_UtX)
    G_each = dat_input$YMMYt - YMU_tauD_UtU_inv_Ut %*% dat_input$MYt - left %*% middle %*% right
    log_det_tauK_I = determinant(1/tau*diag(1/dat_input$delta)+ dat_input$UtU, logarithm=TRUE)$modulus[1] + determinant(tau*diag(dat_input$delta), logarithm=TRUE)$modulus[1]
    Xt_invmiddle_X = dat_input$XtX - dat_input$XtU %*% solve(dat_input$UtU + 1/tau *diag( 1/dat_input$delta)) %*% dat_input$UtX
    log_det_Xt_inv_X = determinant(Xt_invmiddle_X, logarithm=TRUE)$modulus[1]

    sum_det=0
    sum_det=sum_det+(0.5*log_det_tauK_I+0.5*log_det_Xt_inv_X  )*PCnum
    W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
    -(-sum_det -(k*(n-q))/2*log(dat_input$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}

    



#' @title Estimate parameters in SpatialPCA.
#' @description This function estimates parameters through optimization.
#' @param maxiter: Maximum iteration number.
#' @param log_tau_ini: Initial value of log(tau). Because we need tau to be positive, we calculate exp(log(tau)) during iterations.
#' @param dat_input: Object created in the data_prepare_func function.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects
#' \item{par}{estimated log(tau) value}
#' \item{value}{-log likelihood corresponding to the estimated log(tau) value}
#' \item{counts}{number of calls of the SpatialPCA_estimate function}
#' \item{convergence}{An integer code. 0 indicates successful completion, 1 indicates iteration limit maxiter had been reached.}
#' \item{message}{A character string giving any additional information returned by the optimizer}
#' \item{T_m_trend}{Total running time for this function}
#' @export
SpatialPCA_estimate_parameter = function(maxiter=300,log_tau_ini=0,dat_input,PCnum=PCnum){
  suppressMessages(require(RSpectra))
  set.seed(1234)
  param_ini=log_tau_ini
  start_time <- Sys.time()
  m_trend=try(optim(param_ini, SpatialPCA_estimate,dat_input=dat_input,PCnum=PCnum,control = list(maxit = maxiter), lower = -10, upper = 10,method="Brent"),silent=T)
  end_time <- Sys.time()
  T_m_trend = end_time - start_time
  T_m_trend
  m_trend$T_m_trend = T_m_trend

return(m_trend)
}





#' @title Estimate the loading matrix.
#' @description This function estimates loading matrix W in SpatialPCA.
#' @param parameter: The estimated log(tau) value.
#' @param dat_input: The input data object.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects.
#' \item{W}{The estimated loading matrix W.}
#' \item{sigma2_0}{The estimated sigma^2_0 value.}
#' @export
SpatialPCA_estimate_W = function(parameter, dat_input,PCnum=PCnum){
    suppressMessages(require(RSpectra))

    param = parameter
    tau=exp(param[1])
    k = dim(dat_input$Y)[1]
    n = dim(dat_input$Y)[2]
    q=dat_input$q
    
    tauD_UtU_inv = solve(tau*diag(dat_input$delta) + dat_input$UtU)
    YMU_tauD_UtU_inv_Ut = dat_input$YMU %*% tauD_UtU_inv %*% dat_input$Ut
    YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% dat_input$H
    XtU_inv_UtX = dat_input$XtU %*% tauD_UtU_inv %*% dat_input$UtX
    left = dat_input$YMX - YMU_tauD_UtU_inv_UtX
    right = t(left)
    middle = solve(-XtU_inv_UtX)
    G_each = dat_input$YMMYt - YMU_tauD_UtU_inv_Ut %*% dat_input$MYt - left %*% middle %*% right
    W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors

    return_list=list(1:2)
    return_list[[1]]=W_est_here
    return_list[[2]]=(dat_input$tr_YMY+F_funct_sameG(W_est_here,G_each))/(k*(n-q))
    names(return_list) = c("W","sigma2_0")
    return_list

}





#' @title Estimate spatial PCs in SpatialPCA.
#' @description This function estimates spatial PC matrix Z.
#' @param parameter: The estimated log(tau) value.
#' @param dat_input: The input data object.
#' @param estW: The output from function SpatialPCA_estimate_W.
#' @param PCnum: Number of spatial PCs.
#' @param fast: A logic value, TRUE if one wants to use approximation to quickly calculate the spatial PCs, otherwise select FALSE. It is recommended to set fast=TRUE when the sample size is large.
#' @param eigenvecnum: An integer, the number of top eigenvectors and eigenvalues to be used in low-rank approximation in the eigen decomposition step for kernel matrix. 
#' The default is NA, if specified, it is recommended that these top eigen values explain 90% of the variance.
#' @return A list of objects.
#' \item{Z_hat}{The estimated spatial PC matrix Z matrix.} 
#' \item{U}{Eigenvectors used in K.} 
#' \item{delta}{Eigenvalues used in K.} 
#' @export
SpatialPCA_estimate_Z = function(parameter,dat_input,estW,PCnum,fast=FALSE,eigenvecnum=NA){

suppressMessages(require(RSpectra))

    n = dim(dat_input$Y)[2]
    Z_hat = matrix(0, PCnum, n)
    tau = exp(parameter[1])
    W_hat = estW[[1]]
    sigma_2_0_here = estW[[2]]
    sigma_2_0_here = as.numeric(sigma_2_0_here)

if(fast==FALSE){
        U=dat_input$U
        delta=dat_input$delta
}else if(fast==TRUE){
    
    if(is.na(eigenvecnum)==FALSE){
        print(paste0("Low rank approximation!"))
        print(paste0("Using user selected top ",eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
        EIGEN = eigs_sym(dat_input$K, k=eigenvecnum, which = "LM")
        U=EIGEN$vectors
        delta=EIGEN$values

    }else if(n>5000){
        fast_eigen_num = ceiling(n*0.1)
        print(paste0("Low rank approximation!"))
        print("Large sample, using top 10% sample size of eigenvectors and eigenvalues in the Kernel matrix!")
        EIGEN = eigs_sym(dat_input$K, k=fast_eigen_num, which = "LM")
        U=EIGEN$vectors
        delta=EIGEN$values
    }else{
        U=dat_input$U
        delta=dat_input$delta
        ind=length(delta)
        print(paste0("Low rank approximation!"))
        print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
    }
}

    W_hat_t = t(W_hat)
    WtYM = W_hat_t%*% dat_input$YM
    WtYMK = WtYM %*% dat_input$K
    WtYMU = WtYM %*% U
    Ut=t(U)
    UtM = Ut %*% dat_input$M
    UtMK = UtM %*% dat_input$K
    UtMU = UtM %*% U
    middle_inv = solve(1/tau * diag(1/delta) + UtMU)
    Z_hat = tau*WtYMK - tau*WtYMU %*% middle_inv %*% UtMK

    return(list(Z_hat = Z_hat, U =  U, delta = delta))
}



#' @export
F_funct_sameG = function(X,G){ # G is a matrix
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G%*%X[,i]
  }
  -return_val
}


#' @title High-resolution spatial map construction.
#' @description This function predicts spatial PC values at new spatial locations.
#' @param dat_input: The input data object.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects.
#' \item{Z_star}{Predicted Z matrix on new locations.} 
#' \item{Location_star}{Coordinates of new locations.} 
#' @export
high_resolution = function(dat_input,PCnum){

    info = dat_input$location
    K = dat_input$K
    kernelpara = dat_input$kerneltype
    ED = dat_input$ED
    est_log_tau = dat_input$Est_para$par
    est_W = dat_input$Est_W[[1]]
    est_sigma0 = dat_input$Est_W[[2]][1,1]
    est_Z = dat_input$Est_Z
    beta = bandwidth = dat_input$bandwidth
    n=dim(info)[1]
    tau = exp(est_log_tau)

    dis = c()
    for(i in 1:dim(ED)[1]){
        dis[i] = min(ED[i,-i]) # not count the cell it self, only use distance to its nearest cell
    }
    small_distance = median(dis)/4

    num_obs=dim(info)[1]
    Z_hat=matrix(0,PCnum,num_obs)
    num_obs_all = num_obs
    z_star=matrix(0,PCnum,num_obs_all)
    info_new = matrix(0,num_obs_all,2)

    info_new1 = info_new2 = info_new3 = info_new4 = info

    info_new1[,1] = info[,1] - small_distance
    info_new1[,2] = info[,2] + small_distance
    info_new2[,1] = info[,1] + small_distance
    info_new2[,2] = info[,2] + small_distance
    info_new3[,1] = info[,1] - small_distance
    info_new3[,2] = info[,2] - small_distance
    info_new4[,1] = info[,1] + small_distance
    info_new4[,2] = info[,2] - small_distance

    newinfo = rbind(info_new1,info_new2,info_new3,info_new4)
    num_obs_all = dim(newinfo)[1]

    colnames(info) = c("adj_x","adj_y")
    colnames(newinfo) = c("adj_x","adj_y")
    info_all = rbind(info,newinfo)

    ED_all <- as.matrix(dist(scale(info_all[ ,1:2])))
    ED2_all = ED_all^2

    if(kernelpara == "gaussian"){
        K_all = exp(-1.0*ED2_all/as.numeric(beta))
    }else if(kernelpara == "cauchy"){
        K_all = 1.0/(1.0+ED2_all/as.numeric(beta))
    }else if(kernelpara == "quadratic"){
        K_all =  1.0-ED2_all/(ED2_all+as.numeric(beta))
    }

    Sigma_YX = K_all[(num_obs+1):(num_obs+num_obs_all),1:num_obs]
    K_inv = dat_input$U %*% diag(1/dat_input$delta) %*% t(dat_input$U)
    z_star_t= Sigma_YX %*% K_inv %*% t(dat_input$Est_Z$Z_hat)
    z_star = t(z_star_t)
    
    return(list("Z_star"=z_star, "Location_star" = newinfo))
}

