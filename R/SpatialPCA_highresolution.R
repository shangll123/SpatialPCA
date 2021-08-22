########################################################################################################################
# Package: SpatialPCA
# Version: 1.0.1
# Date   : 2021-07-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics 
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu 
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' High-resolution spatial map construction.
#' @param object: SpatialPCA object.
#' @return Returns SpatialPCA object with estimated Spatial PCs.
#' @export

SpatialPCA_highresolution = function(object){

    info = scale(object@location)
    K = object@kernelmat
    ED = 1*as.matrix(dist(info))
    tau = object@tau
    est_W = object@W
    est_sigma0 = object@sigma2_0
    n=dim(info)[1]
    
    #-------------------------------------------------------------------------------#
    # determine small distance to jitter the center cell to surrounding 4 cells		#
    #-------------------------------------------------------------------------------#
    dis = c()
    for(i in 1:dim(ED)[1]){
        dis[i] = min(ED[i,-i]) # not count the cell it self, only use distance to its nearest cell
    }
    small_distance = median(dis)/4

    z_star=matrix(0,object@SpatialPCnum,n)
    info_new = matrix(0,n,2)
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
    # ED2_all <- as.matrix(dist(scale(info_all[ ,1:2]))^2)

    #-------------------------------------------------------------------------------#
    # calculate kernel matrix from new locations									#
    #-------------------------------------------------------------------------------#


	if (object@kerneltype == "gaussian") {
        K_all = exp(-1*as.matrix(dist(info_all)^2)/object@bandwidth)
    }else if (object@kerneltype == "cauchy") {
        K_all = 1/(1 + 1*as.matrix(dist(info_all)^2)/as.numeric(object@bandwidth))
    }else if (object@kerneltype == "quadratic") {
    	ED2=1*as.matrix(dist(info_all)^2)
        K_all = 1 - ED2/(ED2 + as.numeric(object@bandwidth))
    }


    Sigma_YX = K_all[(n+1):(n+num_obs_all),1:n]
    K_inv = object@params$U %*% diag(1/object@params$delta) %*% t(object@params$U)

    z_star_t= Sigma_YX %*% K_inv %*% t(object@SpatialPCs)
    z_star = t(z_star_t)

    rownames(newinfo) = NULL
    object@highPCs = z_star
    object@highPos = newinfo

    return(object)
}









