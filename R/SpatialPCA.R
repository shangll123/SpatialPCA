

#' @title Scale expression values of each gene
#' @description This function scales the expression of each input gene.
#' @param expr: A m gene by n location of gene expression matrix.
#' @return A m gene by n location of gene expression matrix with each row of gene scaled
#' @export
scale_expr = function(expr){
print(paste0("Input expression data: ",dim(expr)[1]," genes on ",dim(expr)[2]," locations."))
  for(i in 1:dim(expr)[1]){
  expr[i,] = scale(expr[i,])
}
  expr
}






#' @title Calculate and prepare necessary matrices, such as projection matrix, distance matrix.
#' @description This function calculates necessary objects before applying SpatialPCA.
#' @param expr: A m gene by n location of gene expression matrix.
#' @param info: A n cell by k dimension of location matrix. n is cell number, k is dimension (k=2 when location is 2D; k=3 when location is 3D)
#' @param covariate: A covariate by n location matrix. If no location related covariate used, the default is "NA". 
#' @return A list of data matrices.
#' \item{YM}{YM matrix}
#' \item{ED}{ED distance matrix}
#' \item{ED2}{ED2 distance square matrix}
#' \item{tr_YMY}{tr_YMY scalar}
#' \item{M}{M matrix}
#' \item{H}{H matrix}
#' \item{n}{n matrix}
#' \item{Y}{expression matrix}
#' \item{q}{number of covariates, if no input covariate is avaliable then q=1 because we added an intercept in the function.}
#' @export
data_prepare_func = function(expr,info,covariate=NA){

  if (ncol(expr) != nrow(info)) {
    stop("ERROR - expression samples size doesn't agree with the number of locations")
  }

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
res<-list("YM"=YM, "ED"=ED, "ED2"=ED2,  "tr_YMY"=tr_YMY, "M"=M, "H"=H, "n"=n,"Y"=expr,"q"=q)
return(res)

}





#' @title Select bandwidth in Gaussian kernel
#' @description This function selects bandwidth in Gaussian kernel.
#' @param expr: A m gene by n location of gene expression matrix.
#' @param info: A n cell by k dimension of location matrix. n is cell number, k is dimension 
#' @param method: Method used in bandwidth selection, "SJ" for Sheather & Jones (1991) method, 
#' "Silverman" for Silverman's ‘rule of thumb’ method (1986), "Scott" for Scott (1992) method
#' @return A scalar of calculated bandwidth.
#' @export

bandwidth_select=function (expr, info, method) 
{
    N = dim(info)[1]
    if (method == "SJ") {
        
              bw_SJ = c()
        for (i in 1:dim(expr)[1]) {
            tryCatch({ 
              print(i)
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
#' @description This function selects bandwidth in Gaussian kernel.
#' @param kernelpara: An integer that represents different kernels: "gaussian", "cauchy", or "quadratic".
#' @param ED2: An n by n distance square matrix.
#' @param beta: A scalar of bandwidth value.
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






#' @title Calculate -log likelihood value in parameter estimation.
#' @description This function calculate -log likelihood value.
#' @param param_ini: Input log(tau) value. 
#' @param dat_input: Object obtained from data_prepare_func function
#' @param PCnum: Number of spatial PCs.
#' @return A numeric value of -log likelihood.
#' @export
SpatialPCA_estimate = function(param_ini,dat_input,PCnum){
   param = param_ini
    tau=exp(param[1])
    print(param[1])
    k = dim(dat_input$Y)[1]
    n = dim(dat_input$Y)[2]
    sum_det=0
    q=dat_input$q 
    Sigma=tau*K
    Sigma_tilde=Sigma+diag(dat_input$n)
    L=t(chol(Sigma_tilde)) 
    L_H=t(chol(t(dat_input$H) %*% solve(Sigma_tilde) %*% dat_input$H))
    K_inv = solve(K)
    G_each = dat_input$Y %*% dat_input$M %*% solve(dat_input$M +1/tau*K_inv)  %*% dat_input$M%*%t(dat_input$Y)
    for(i_d in 1:PCnum){
    sum_det=sum_det+(sum(log(diag(L)))+sum(log(diag(L_H))))
  }
  W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
  -(-sum_det -(k*(n-q))/2*log(dat_input$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}






#' @title Estimate parameters in SpatialPCA.
#' @description This function estimates parameter tau with optim function.
#' @param maxiter: Maximum iteration number.
#' @param log_tau_ini: Initial value of log(tau). Because we want to make sure tau is greater than 0 during iterations.
#' @param dat_input: Object obtained from data_prepare_func function.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects
#' \item{par}{estimated log(tau) value}
#' \item{value}{-log likelihood corresponding to estimated log tau value}
#' \item{counts}{number of calls to the SpatialPCA_estimate function}
#' \item{convergence}{An integer code. 0 indicates successful completion, 1 indicates iteration limit maxiter had been reached.}
#' \item{message}{A character string giving any additional information returned by the optimizer}
#' \item{T_m_trend}{Total running time for this function}
#' @export
SpatialPCA_estimate_parameter = function(maxiter=300,log_tau_ini=0,dat_input,PCnum=20){
  suppressMessages(require(RSpectra))
  set.seed(1234)
  param_ini=log_tau_ini
  start_time <- Sys.time()
  m_trend=try(optim(param_ini, SpatialPCA_estimate,dat_input=dat_input,PCnum=PCnum,control = list(maxit = maxiter), lower = -3, upper = 3,method="Brent"),silent=T)
  end_time <- Sys.time()
  T_m_trend = end_time - start_time
  T_m_trend
  m_trend$T_m_trend = T_m_trend

return(m_trend)
}





#' @title Estimate W matrix in SpatialPCA.
#' @description This function estimates loading matrix W.
#' @param parameter: Estimated log(tau) value.
#' @param dat_input: Object obtained from data_prepare_func function.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects, the first one is estimated W matrix under given parameter tau; the second one is estimated sigma_0 value.
#' @export
SpatialPCA_estimate_W = function(parameter, dat_input,PCnum=20){
   param = parameter
    tau=exp(param[1])
    Sigma=tau*K
    Y = dat_input$Y
    n=dim(dat_input$Y)[2]
    q=dat_input$q
    k=dim(dat_input$Y)[1]
    Sigma_tilde=Sigma+diag(n)
    L=t(chol(Sigma_tilde))
    L_H=t(chol(t(dat_input$H)%*%solve(Sigma_tilde)%*% dat_input$H))
    sum_det=0
    K_inv = solve(K)
    G_each = dat_input$Y %*% dat_input$M %*% solve(dat_input$M +1/tau*K_inv)  %*% dat_input$M%*%t(dat_input$Y)
    for(i_d in 1:PCnum){
    sum_det=sum_det+(sum(log(diag(L)))+sum(log(diag(L_H))))
  }
  W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
  return_list=list(1:2)
  return_list[[1]]=W_est_here
  return_list[[2]]=(dat_input$tr_YMY+F_funct_sameG(W_est_here,G_each))/(k*(n-q))
  return_list

}







#' @title Estimate Z matrix in SpatialPCA.
#' @description This function estimates spatial PC matrix Z.
#' @param parameter: Estimated log(tau) value.
#' @param dat_input: Object obtained from data_prepare_func function.
#' @param estW: Estimated W matrix.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects
#' \item{Z_hat}{estimated Z matrix} 
#' \item{mat_inv}{mat_inv matrix, will be used in high-resolution prediction} 
#' \item{YM_mat_inv}{YM_mat_inv matrix, will be used in high-resolution prediction} 
#' @export
SpatialPCA_estimate_Z = function(parameter,dat_input,estW,PCnum=20){
n = dim(dat_input$Y)[2]
    Z_hat = matrix(0, PCnum, n)
    tau = exp(parameter[1])
    W_hat = estW[[1]]
    sigma_2_0_here = estW[[2]]
    sigma_2_0_here = as.numeric(sigma_2_0_here)

    mat_inv = solve( diag(n)+tau *  K %*% dat_input$M  )
    YM_mat_inv = dat_input$YM %*% mat_inv
    for (i_d in 1:PCnum) {
        print(i_d)
        middle_part = t(W_hat[, i_d]) %*% YM_mat_inv
        middle_part = as.matrix(middle_part)
        Z_hat[i_d, ] = tau * middle_part %*%   K
    }
    return(list(Z_hat = Z_hat, YM_mat_inv =  1/sigma_2_0_here*YM_mat_inv))
}


#' @export
F_funct_sameG = function(X,G){ # G is a matrix
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G%*%X[,i]
  }
  -return_val
}


#' @title Calculate -log likelihood value in parameter estimation in large sample size.
#' @description This function calculate -log likelihood value.
#' @param param_ini: Input log(tau) value. 
#' @param dat_input: Object obtained from data_prepare_func function
#' @param PCnum: Number of spatial PCs.
#' @return A numeric value of -log likelihood.
#' @export
SpatialPCA_estimate_paras_largedata = function(param_ini,dat_input,PCnum=20){

  suppressMessages(require(RSpectra))
    set.seed(1234)
    param = param_ini
    print(param)
    tau=exp(param[1])
    k = dim(dat_input$Y)[1]
    n = dim(dat_input$Y)[2]
    q=dat_input$q

# # Q = (\tau *K+I)^-1
    UtU = t(dat_input$U) %*% dat_input$U
    Q_middle_inv = solve(UtU + 1/tau *diag( 1/dat_input$delta)) # d by d
    Q = diag(n)-dat_input$U %*% Q_middle_inv %*% t(dat_input$U ) # n by n
    XtQ = t(dat_input$H) %*% Q # q by n
    XtQX = XtQ %*% dat_input$H
# P=(I+\tau^-1 K^-1)^-1
    P_middle_inv = solve(UtU + tau *diag( dat_input$delta)) # d by d
    P = diag(n)-dat_input$U %*% P_middle_inv %*% t(dat_input$U ) # n by n
    XtP = t(dat_input$H) %*% P # q by n
    XtPX = XtP %*% dat_input$H 
    mat_middle_inv = solve(-t(dat_input$H) %*% dat_input$H +XtPX ) # q by q
    W_middle_inv = P - t(XtP) %*% mat_middle_inv %*% XtP # n by n
# W_middle_inv = (M+\tau^-1 K^-1)^-1
    G_each = as.matrix(dat_input$YM%*% W_middle_inv %*% dat_input$YMt ) # m by m
# det(\tau *K +I) = det(\tau^-1 *D^-1 + UtU)*det(\tau *D)
    det_tauK_I = det(1/tau*diag(1/dat_input$delta)+ UtU) * det(tau*diag(dat_input$delta))
    sum_det=0
    sum_det=sum_det+(0.5*log(det_tauK_I)+0.5*log(det(XtQX))  )*PCnum
    W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
    -(-sum_det -(k*(n-q))/2*log(dat_input$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}








#' @title Estimate parameters in SpatialPCA when sample size is large.
#' @description This function estimates parameters through optim function.
#' @param maxiter: Maximum iteration number.
#' @param log_tau_ini: Initial value of log(tau).
#' @param dat_input: Object obtained from data_prepare_func function
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects
#' \item{par}{estimated tau value}
#' \item{value}{-log likelihood corresponding to estimated tau value}
#' \item{counts}{number of calls to the SpatialPCA_estimate function}
#' \item{convergence}{An integer code. 0 indicates successful completion, 1 indicates iteration limit maxiter had been reached.}
#' \item{message}{A character string giving any additional information returned by the optimizer}
#' \item{T_m_trend}{Total running time for this function}
#' @export
SpatialPCA_estimate_parameter_largedata = function(maxiter=300,log_tau_ini=0,dat_input,PCnum=20){

  suppressMessages(require(RSpectra))
  set.seed(1234)
  param_ini=log_tau_ini
  start_time <- Sys.time()
  m_trend=try(optim(param_ini, SpatialPCA_estimate_paras_largedata,dat_input=dat_input,PCnum=PCnum,control = list(maxit = maxiter), lower = -10, upper = 10,method="Brent"),silent=T)
  end_time <- Sys.time()
  T_m_trend = end_time - start_time
  T_m_trend
  m_trend$T_m_trend = T_m_trend

return(m_trend)
}






#' @title Estimate W matrix in SpatialPCA when sample size is large.
#' @description This function estimates loading matrix W.
#' @param parameter: Estimated log(tau) value.
#' @param dat_input: Object obtained from data_prepare_func function
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects, the first one is estimated W matrix under given parameter tau; the second one is estimated sigma_0 value.
#' @export
SpatialPCA_estimate_W_largedata = function(parameter,dat_input,PCnum=20){

 suppressMessages(require(RSpectra))

    param = parameter
    tau=exp(param[1])

    k = dim(dat_input$Y)[1]
    n = dim(dat_input$Y)[2]
    q=dat_input$q
    UtU = t(dat_input$U) %*% dat_input$U
# P=(I+\tau^-1 K^-1)^-1
    P_middle_inv = solve(UtU + tau *diag( dat_input$delta)) # d by d
    P = diag(n)-dat_input$U %*% P_middle_inv %*% t(dat_input$U ) # n by n
    XtP = t(dat_input$H) %*% P # q by n
    XtPX = XtP %*% dat_input$H
    mat_middle_inv = solve(-t(dat_input$H) %*% dat_input$H +XtPX ) # q by q
    W_middle_inv = P - t(XtP) %*% mat_middle_inv %*% XtP # n by n
# W_middle_inv = (M+\tau^-1 K^-1)^-1
    G_each = as.matrix(dat_input$YM%*% W_middle_inv %*% dat_input$YMt )
    W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
    return_list=list(1:2)
    return_list[[1]]=W_est_here
    return_list[[2]]=(dat_input$tr_YMY+F_funct_sameG(W_est_here,G_each))/(k*(n-q))
    return_list
}





#' @title Estimate Z matrix in SpatialPCA when sample size is large.
#' @description This function estimates spatial PC matrix Z.
#' @param parameter: Estimated log(tau) value.
#' @param dat_input: Object obtained from data_prepare_func function.
#' @param estW: Estimated W matrix.
#' @param PCnum: Number of spatial PCs.
#' @return A list of objects
#' \item{Z_hat}{estimated Z matrix} 
#' \item{mat_inv}{mat_inv matrix, will be used in high-resolution prediction} 
#' \item{YM_mat_inv}{YM_mat_inv matrix, will be used in high-resolution prediction} 
#' @export
SpatialPCA_estimate_Z_largedata = function(parameter,dat_input,estW,PCnum=20 ){
n = dim(dat_input$Y)[2]
    Z_hat = matrix(0, PCnum, n)
    tau = exp(parameter[1])
    W_hat = estW[[1]]
    sigma_2_0_here = estW[[2]]
    sigma_2_0_here = as.numeric(sigma_2_0_here)

    mat_inv = solve(tau * K %*% dat_input$M +  diag(n))
    YM_mat_inv = dat_input$YM %*% mat_inv
    for (i_d in 1:PCnum) {
        print(i_d)
        middle_part = t(W_hat[, i_d]) %*% YM_mat_inv
        middle_part = as.matrix(middle_part)
        Z_hat[i_d, ] = tau * middle_part %*%   K
    }
    return(list(Z_hat = Z_hat, mat_inv = 1/sigma_2_0_here*mat_inv, YM_mat_inv = YM_mat_inv))
}







#' @title High-resolution spatial map construction.
#' @description This function predicts spatial PC values at new spatial locations.
#' @param info: A n by k matrix, current location matrix.
#' @param K: A n by n matrix, current kernel matrix.
#' @param kernelpara: An integer that represents different kernels: 1 - Gaussian kernel, 2 - cauchy kernel, 3 - quadratic kernel. Here it should be consistent with the kernel that used in modeling original spatial PCs.
#' @param ED: A n by n Distance matrix.
#' @param est_log_tau: Estimated log tau value.
#' @param est_W: Estimated loading matrix W. 
#' @param est_sigma0: Estimated sigma_0 value.
#' @param est_Z: A object from SpatialPCA_estimate_Z function.
#' @param PCnum: Number of predicted spatial PCs.
#' @return A list of objects
#' \item{Z_star}{Predicted Z matrix on new locations} 
#' \item{Location_star}{Coordinates of new locations} 
#' @export
high_resolution = function(info, K, kernelpara, ED,est_log_tau, est_W, est_sigma0, est_Z,PCnum){
#high_resolution = function(info, K, result,est_Z){
n=dim(info)[1]
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


beta=bandwidth
tau = exp(est_log_tau)

if(kernelpara == "gaussian"){
    K_all = exp(-1.0*ED2_all/as.numeric(beta))
}else if(kernelpara == "cauchy"){
    K_all = 1.0/(1.0+ED2_all/as.numeric(beta))
}else if(kernelpara == "quadratic"){
    K_all =  1.0-ED2_all/(ED2_all+as.numeric(beta))
}


Z_hat=matrix(0,PCnum,n)
W_hat = est_W
sigma_2 = tau*est_sigma0
Sigma=sigma_2*K

Sigma_all = sigma_2*K_all
Sigma_YX = Sigma_all[(num_obs+1):(num_obs+num_obs_all),1:num_obs]

z_star=matrix(0,PCnum,num_obs_all)

      for(i_d in 1:PCnum){
        print(i_d)
        middle_part=t(W_hat[,i_d])%*%est_Z$YM_mat_inv
        middle_part = as.matrix(middle_part)
        z_star[i_d,]=Sigma_YX %*% t(middle_part)
    }

return(list("Z_star"=z_star, "Location_star" = newinfo))
}
