########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' High-resolution spatial map construction.
#' @param object SpatialPCA object.
#' @param platform 
#' "ST": the 10X ST platform, impute 9 subspots in a square shape for each spot; 
#' "Visium": the 10X Visium platform, impute 6 subspots in a hexagonal shape for each spot; 
#  "Other": impute data from any platform, jitter each existing location to surrounding 4 locations as new locations.
#  if platform is selected as NULL: the user may input newlocation matrix with coordinates on the locations they want to impute.
#' @param newlocation A n* by d location matrix, n* is number of new locations, d is dimension of locations.
#' Users can optionally provide new locations at the original location scale.
#' @return Returns SpatialPCA object with estimated Spatial PCs on new locations.
#' @export

SpatialPCA_highresolution = function(object,platform="ST",newlocation=NULL){

    info = scale(object@location)
    #K = object@kernelmat
    ED = 1*as.matrix(dist(info))
    tau = object@tau
    est_W = object@W
    est_sigma0 = object@sigma2_0
    n=dim(info)[1]

if(is.null(newlocation)){
        if(platform=="ST"){

            newinfo = .make_subspot_coldata(object@location, platform="ST")[,5:6]
            newinfo=as.matrix(newinfo)
            colnames(info) = c("adj_x","adj_y")
            colnames(newinfo) = c("adj_x","adj_y")
            info = object@location
            info_all = scale(rbind(info,newinfo))
            num_obs_all = dim(newinfo)[1]

        }else if(platform=="Visium"){

            newinfo = .make_subspot_coldata(object@location, platform="Visium")[,5:6]
            newinfo=as.matrix(newinfo)
            info = object@location
            colnames(info) = c("adj_x","adj_y")
            colnames(newinfo) = c("adj_x","adj_y")
            info_all = scale(rbind(info,newinfo))
            num_obs_all = dim(newinfo)[1]

        }else if(platform=="Other"){

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

        }
}else if(!is.null(newlocation)){

    newinfo = as.matrix(newlocation)
    num_obs_all = dim(newinfo)[1]
    info = object@location
    colnames(info) = c("adj_x","adj_y")
    colnames(newinfo) = c("adj_x","adj_y")
    info_all = scale(rbind(info,newinfo))

}

    #-------------------------------------------------------------------------------#
    # calculate kernel matrix from new locations                                    #
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




#' High-resolution gene expression prediction.
#' @param object SpatialPCA object with high resolution predicted spatial PCs.
#' Users can optionally provide new locations at the original location scale.
#' @return Returns SpatialPCA object with predicted gene expression on new locations.
#' We can predict normalized gene expression for the genes in the object@normalized_expr matrix.
#' @export
SpatialPCA_expr_pred = function(object){
    if(!exists(object@highPCs)){
      print("Please first use SpatialPCA_highresolution function to predict high resolution spatial PCs.")
      return(object)
    }else{
    object@expr_pred =  object@W %*% object@highPCs
    return(object)
    }
}


# follow BayesSpace to make subspots
.make_subspot_offsets <- function(n_subspots_per) {
    if (n_subspots_per == 6) {
        rbind(expand.grid(c(1/3, -1/3), c(1/3,-1/3)), expand.grid(c(2/3, -2/3), 0))
    # } else if (n_subspots_per == 7) {
    #     rbind(expand.grid(c(1/3, -1/3), c(1/3, -1/3)), expand.grid(c(2/3, -2/3, 0), 0))
    } else if (n_subspots_per == 9) {
        rbind(expand.grid(c(1/3, -1/3, 0), c(1/3, -1/3, 0)))
    } else {
        stop("Only 6 and 9 subspots currently supported.")
    }
}

# follow BayesSpace to make subspots
.make_subspot_coldata <- function(positions, platform="ST") {

    require(assertthat)
    n_subspots_per <- ifelse(platform == "Visium", 6, 9)

    positions <- as.data.frame(positions)
    colnames(positions)=c("row","col")
    #colnames(cdata) <- c("imagecol", "imagerow")
    
    n_spots <- nrow(positions)
    n_subspots <- n_spots*n_subspots_per
    cdata = data.frame("row"=rep(NA,n_subspots),"col"=rep(NA,n_subspots))
    assertthat::assert_that(nrow(cdata) == n_spots * n_subspots_per)
    
    ## Index of parent spot is (subspot % n_spots)
    idxs <- seq_len(n_subspots)
    spot_idxs <- ((idxs - 1) %% n_spots) + 1
    subspot_idxs <- rep(seq_len(n_subspots_per), each=n_spots)
    cdata$spot.idx <- spot_idxs
    cdata$subspot.idx <- subspot_idxs
    rownames(cdata) <- paste0("subspot_", spot_idxs, ".", subspot_idxs)
    
    offsets <- .make_subspot_offsets(n_subspots_per)
    cdata$spot.row <- rep(positions$row, n_subspots_per)
    cdata$spot.col <- rep(positions$col, n_subspots_per)
    cdata$col <- cdata$spot.col + rep(offsets[, 1], each=n_spots)
    cdata$row <- cdata$spot.row + rep(offsets[, 2], each=n_spots)

    cols <- c("spot.idx", "subspot.idx", "spot.row", "spot.col", "row", "col")
    cdata[, cols]
}












