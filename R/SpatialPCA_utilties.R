########################################################################################################################
# Package: SpatialPCA
# Version: 1.1.0
# Date   : 2021-10-27
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics
# Authors: L. Shang and X. Zhou
# Contacts: shanglu@umich.edu
#          University of Michigan, Department of Biostatistics
####################################################################################################


#' Obtain PCA low dimensional components for dimension reduction methods comparison.
#' @param expr Normalized gene expression g by n matrix. g is gene number, n is sample size.
#' @param PCnum Number of PCs.
#' @return A d by n matrix of low dimensional components from PCA. d is number of low dimensional components, n is sample size.
#' @export
get_PCA = function(expr,PCnum){

  n=dim(expr)[2]
  k = dim(expr)[1]
  output_sub_mean=matrix(0,k,n)
  for(i_k in 1:k){
    output_sub_mean[i_k,]=expr[i_k,]-mean(expr[i_k,])
  }
  svd_output_sub_mean=svd(output_sub_mean)
  A_ini=svd_output_sub_mean$u[,1:PCnum]
  Z_pca = t(A_ini) %*% output_sub_mean
  return(Z_pca)
}

#' Obtain NMF low dimensional components for dimension reduction methods comparison.
#' @param count Count expression g by n matrix. g is gene number, n is sample size.
#' @param PCnum Number of PCs.
#' @return A d by n matrix of low dimensional components from NMF. d is number of low dimensional components, which is same as number of Spatial PCs in SpatialPCA. n is sample size.
#' @export
get_NMF = function(count, PCnum){
  #suppressMessages(require(scater))
  #suppressMessages(require(NMF))
  expr = log(count+1) # non negative
  res <- scater::calculateNMF(expr, ncomponents = PCnum)
  Z_NMF = t(res)
  return(Z_NMF)
}



#' @title Obtain clustering cluster labels through louvain method.
#' @description This function performs louvain clustering on input low dimensional components.
#' @param clusternum The desired number of clusters the user wants to obtain.
#' @param latent_dat A d by n matrix of low dimensional components, d is number of PCs, n is number of spots.
#' @param knearest An integers, number of nearest neighbors for KNN graph construction in louvain clustering.
#' @return The cluster labels.
#'
#' @export
louvain_clustering = function(clusternum, latent_dat,knearest=100){
  set.seed(1234)
  # suppressMessages(require(FNN))
  # suppressMessages(require(igraph))
  # suppressMessages(require(bluster))

  PCvalues = latent_dat
  info.spatial = as.data.frame(t(PCvalues))
  colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))
  knn.norm = FNN::get.knn(as.matrix(t(PCvalues)), k = knearest)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
  k=knearest), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = igraph::simplify(nw.norm)
  lc.norm = igraph::cluster_louvain(nw.norm)
  merged <- bluster::mergeCommunities(nw.norm, lc.norm$membership, number=clusternum)
  clusterlabel = as.character(as.integer(as.factor(paste0("cluster",merged))))
  return("cluster_label"=clusterlabel)
}


#' @title Obtain clustering cluster labels through walktrap method.
#' @description This function performs walktrap clustering on input low dimensional components.
#' @param clusternum The desired number of clusters the user wants to obtain.
#' @param latent_dat A d by n matrix of low dimensional components, d is number of PCs, n is number of spots.
#' @param knearest An integers, number of nearest neighbors for SNN graph construction in walktrap clustering.
#' @return The cluster labels.
#'
#' @export
walktrap_clustering = function(clusternum, latent_dat,knearest=100){
  set.seed(1234)
  # suppressMessages(require(bluster))
  # suppressMessages(require(igraph))

  PCvalues = latent_dat
  g <- bluster::makeSNNGraph(as.matrix(t(PCvalues)),k = knearest)
  #clusters <- igraph::cluster_fast_greedy(g)$membership
  g_walk <- igraph::cluster_walktrap(g)
  cluster_label_new = as.character(igraph::cut_at(g_walk, no=clusternum))

  return("cluster_label"=cluster_label_new)
}




#' @title Refine clustering results for 10x ST or Visium data.
#' @description This function refines spatial clustering of ST or Visium data.
#' @param clusterlabels The cluster label obtained (e.g. from louvain method or walktrap method).
#' @param location A n by 2 location matrix of spots.
#' @param shape Select shape='hexagon' for Visium data, 'square' for ST data.
#' @return The refined cluster labels.
#' @export
refine_cluster_10x = function(clusterlabels, location, shape="square"){

  dis_df = as.matrix(dist(location))
  if(shape=="square"){
    num_obs = 4
  }else if(shape == "hexagon"){
    num_obs = 6
  }else{
    print("Select shape='hexagon' for Visium data, 'square' for ST data.")
  }
  refined_pred = clusterlabels
  for(i in 1:length(clusterlabels)){
    nearby_spot_ind = order(dis_df[i,])[1:(num_obs+1)]
    labels_nearby_spot_ind = refined_pred[nearby_spot_ind]  # use updated cluster
    spot_of_interest = refined_pred[i]
    labels_table = table(labels_nearby_spot_ind)
    if( labels_table[spot_of_interest]<num_obs/2 & max(labels_table)>num_obs/2){
      refined_pred[i] = names(labels_table)[which.max(labels_table)]
    }else{
      refined_pred[i] = spot_of_interest
    }

  }

  return(refined_pred)
}




#' @title Visualize PCs on their locations.
#' @description This function visualizes the low dimensional component values.
#' @param location A n cell by k dimension of location matrix. n is cell number, k=2 if the spots are on 2D space.
#' @param PCs A d by n matrix of low dimensional components.
#' @param textmethod A text string of the name of method used to extract latent factors, e.g. "SpatialPCA" or "PCA". It will be shown as the title of the figures.
#' @param pointsize The point size of each location for visualization.
#' @param textsize The text size in the figure legend.
#' @return A list of ggplot objects for factor value plots.
#'
#' @import ggplot2
#'
#' @export
plot_factor_value = function(location, PCs,textmethod,pointsize=2,textsize=15){

  location = as.data.frame(location)
  PCnum=dim(PCs)[1]
  p = list()
  for(k in 1:PCnum){
      locc1 = location[,1]
      locc2 = location[,2]
      PC_value = PCs[k,]
      datt = data.frame(PC_value, locc1, locc2)
      p[[k]] = ggplot(datt, aes(x = locc1, y = locc2, color = PC_value)) +
        geom_point(size=pointsize, alpha = 1) +
        scale_color_gradientn(colours = c("#4E84C4", "#FFDB6D")) +
        ggtitle(paste0(textmethod," PC ",k))+
        theme_void()+
        theme(plot.title = element_text(size = textsize),
              text = element_text(size = textsize),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 22) ,
              legend.position = "bottom")

}
  return(p)
}




#' @title Visualize cluster labels on locations.
#' @description This function visualizes cluster labels on locations.
#' @param location A n by k matrix of spot locations.
#' @param clusterlabel A vector of cluster labels for spots.
#' @param pointsize An integer, the point size of each spot.
#' @param textsize An integer, the text size in the legend.
#' @param title_in A character string, the title you want to display at the top of the figure.
#' @param color_in A vector of colors for each cluster.
#' @param legend A character string, the position of the figure legend. Select from "top", "bottom","left" or "right".
#' @return A ggplot object.
#' @export
plot_cluster = function(location, clusterlabel, pointsize=3,text_size=15 ,title_in,color_in,legend="none"){
  cluster = clusterlabel
  loc_x=location[,1]
  loc_y=location[,2]
  datt = data.frame(cluster, loc_x, loc_y)
  p = ggplot(datt, aes(x = location[,1], y = location[,2], color = cluster)) +
        geom_point( alpha = 1,size=pointsize) +
        scale_color_manual(values = color_in)+
        ggtitle(paste0(title_in))+
        theme_void()+
        theme(plot.title = element_text(size = text_size,  face = "bold"),
              text = element_text(size = text_size),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 15) ,
              legend.position =legend)
  p
}





#' @title Visualize RGB plot from tSNE.
#' @description We summarized the inferred low dimensional components into three tSNE components and visualized the three resulting components with red/green/blue (RGB) colors in the RGB plot.
#' @param location A n by k location matrix. n is spot number.
#' @param latent_dat A d by n matrix of low dimensional components.
#' @param pointsize The point size of each spot.
#' @param textsize The text size in the legend.
#' @return A list.
#' \item{RGB}{A data frame with five columns: x coordinate, y coordinate, R, G, and B color index}
#' \item{figure}{A ggplot object for RGB plot from tSNE}
#'
#' @import Rtsne
#'
#' @export
plot_RGB_tSNE=function(location, latent_dat,pointsize=2,textsize=15){

  # suppressMessages(require(Rtsne))

  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")

  PCvalues = latent_dat

  tsne <- Rtsne(t(PCvalues),dims=3,check_duplicates = FALSE)
  r = (tsne$Y[,1]-min(tsne$Y[,1]))/(max(tsne$Y[,1])-min(tsne$Y[,1]))
  g = (tsne$Y[,2]-min(tsne$Y[,2]))/(max(tsne$Y[,2])-min(tsne$Y[,2]))
  b = (tsne$Y[,3]-min(tsne$Y[,3]))/(max(tsne$Y[,3])-min(tsne$Y[,3]))
  x =  info$sdimx
  y =  info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
      geom_point(size=pointsize) +
      scale_color_identity()+
      ggtitle(paste0("RGB tSNE"))+
      theme_void()+
      theme(plot.title = element_text(size = textsize),
              text = element_text(size = textsize),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 22) ,
              legend.position = "bottom")

  return(list("RGB"=dat,"figure"=p1))
}



#' @title Visualize RGB plot from UMAP.
#' @description We summarized the inferred low dimensional components into three UMAP. components and visualized the three resulting components with red/green/blue (RGB) colors in the RGB plot.
#' @param location A n by k location matrix. n is spot number.
#' @param latent_dat A d by n matrix of low dimensional components.
#' @param pointsize The point size of each spot.
#' @param textsize The text size in the legend.
#' @return A list.
#' \item{RGB}{A data frame with five columns: x coordinate, y coordinate, R, G, and B color index}
#' \item{figure}{A ggplot object for RGB plot from UMAP.}
#'
#' @import umap
#'
#' @export
plot_RGB_UMAP=function(location, latent_dat,pointsize=2,textsize=15){

  # suppressMessages(require(umap))

  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")

  PCvalues = latent_dat

  umap <- umap(t(PCvalues),n_components = 3)
  r = (umap$layout[,1]-min(umap$layout[,1]))/(max(umap$layout[,1])-min(umap$layout[,1]))
  g = (umap$layout[,2]-min(umap$layout[,2]))/(max(umap$layout[,2])-min(umap$layout[,2]))
  b = (umap$layout[,3]-min(umap$layout[,3]))/(max(umap$layout[,3])-min(umap$layout[,3]))
  x =  info$sdimx
  y = info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
      geom_point(size=pointsize) +
      scale_color_identity()+
      ggtitle(paste0("RGB UMAP"))+
      theme_void()+
      theme(plot.title = element_text(size = textsize),
              text = element_text(size = textsize),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 22) ,
              legend.position = "bottom")

  return(list("RGB"=dat,"figure"=p1))
}



#' @title Visualize pseudotimes on locations.
#' @description This function visualizes pseudotimes on locations.
#' @param pseudotime A length n vector of pseudotime inferred from Slingshot.
#' @param location A n by 2 data frame of spot locations.
#' @param clusterlabels A vector of integers, the cluster labels for each spot
#' @param gridnum Number of grids that evenly segment the whole tissue section.
#' @param color_in A vector of character strings representing colors for each cluster.
#' @param pointsize An integer, the point size of each spot
#' @param arrowlength An integer, the length of arrows inside a grid between one spot with smallest pseudotime and largest pseudotime.
#' @param arrowsize An integer, the size of arrows inside a grid between one spot with smallest pseudotime and largest pseudotime.
#' @param textsize An integer, the size of text in the figure.
#' @return A ggplot object.
#' \item{Pseudotime}{A ggplot object visualizing pseudotime on locations.}
#' \item{Arrowplot1}{A ggplot object for arrows pointing from smallest pseudotime and largest pseudotime in each grid.}
#' \item{Arrowplot2}{A ggplot object for arrows pointing from largest pseudotime and smallest pseudotime in each grid.}
#' \item{Arrowoverlay1}{A ggplot object for arrows pointing from smallest pseudotime and largest pseudotime in each grid, overlayed on clustering plot.}
#' \item{Arrowoverlay2}{A ggplot object for arrows pointing from largest pseudotime and smallest pseudotime in each grid, overlayed on clustering plot.}
#' @export
plot_trajectory = function(pseudotime, location,clusterlabels,gridnum,color_in,pointsize=5 ,arrowlength=0.2,arrowsize=1,textsize=22){

  pseudotime_use=pseudotime
  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")
  grids = gridnum

  min_x = min(info$sdimx)
  min_y = min(info$sdimy)
  max_x = max(info$sdimx)
  max_y = max(info$sdimy)

  x_anchor = c()
  for(x_i in 1:(grids+1)){
      space_x = (max_x - min_x)/grids
      x_anchor[x_i] = min_x+(x_i-1)*space_x
  }
  y_anchor = c()
  for(y_i in 1:(grids+1)){
      space_y = (max_y - min_y)/grids
      y_anchor[y_i] = min_y+(y_i-1)*space_y
  }

# label square by num_x, num_y
count = 0
squares = list()
direction_pseudotime_point = list()
start_x_dat = c()
start_y_dat = c()
end_x_dat = c()
end_y_dat = c()
for(num_x in 1:grids){
  for(num_y in 1:grids){

    filter_x = which(info$sdimx >= x_anchor[num_x] & info$sdimx <= x_anchor[num_x+1])
    filter_y = which(info$sdimy >= y_anchor[num_y] & info$sdimy <= y_anchor[num_y+1])
    # find points in each grid
    points_in_grid = intersect(filter_x, filter_y)


    # find min pseudotime and max pseudotime in each grid
    if(length(points_in_grid)>1 & sum(which.min(pseudotime_use[points_in_grid]))>0){
      count = count + 1
      squares[[count]]= intersect(filter_x, filter_y)
      direction_pseudotime_point[[count]] = list()
      direction_pseudotime_point[[count]]$min_point = info[squares[[count]][which.min(pseudotime_use[squares[[count]]])],]
      direction_pseudotime_point[[count]]$max_point = info[squares[[count]][which.max(pseudotime_use[squares[[count]]])],]
      start_x_dat[count] = unlist(direction_pseudotime_point[[count]]$min_point$sdimx)
      start_y_dat[count] = unlist(direction_pseudotime_point[[count]]$min_point$sdimy)
      end_x_dat[count] = unlist(direction_pseudotime_point[[count]]$max_point$sdimx)
      end_y_dat[count] = unlist(direction_pseudotime_point[[count]]$max_point$sdimy)
    }
  }
}


loc1 = info$sdimx
loc2 = info$sdimy

time = pseudotime_use+0.01
datt = data.frame(time, loc1, loc2)
datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
p01=ggplot(datt, aes(x = loc1, y = loc2, color = time)) +
    geom_point( alpha = 1,size=pointsize) +
    scale_color_gradientn(colours = c("red", "green")) +
        theme_void()+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
              text = element_text(size = textsize),
              legend.position = "bottom")
p02= ggplot()+
  geom_segment(aes(x = start_x_dat, y = start_y_dat, xend = end_x_dat, yend = end_y_dat,colour = "black"),
    arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,color="black",data = datt2) +
        theme_void()+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
              text = element_text(size = textsize),
              legend.position = "bottom")

p03= ggplot()+
  geom_segment(aes(x = end_x_dat, y = end_y_dat, xend = start_x_dat, yend = start_y_dat,colour = "black"),
    arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,color="black",data = datt2) +
        theme_void()+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
              text = element_text(size = textsize),
              legend.position = "bottom")

time = pseudotime_use+0.1
datt1 = data.frame(time, loc1, loc2)
p1 = ggplot(datt1, aes(x = loc1, y = loc2)) +
    geom_point( alpha =1,size=pointsize,aes(color=clusterlabels)) +
        theme_void()+
        scale_colour_manual(values=color_in)+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
              text = element_text(size = textsize),
              legend.position = "bottom")
datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
p2= geom_segment(aes(x = start_x_dat, y = start_y_dat, xend = end_x_dat, yend = end_y_dat,colour = "segment"),
    arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,arrow.fill="black",data = datt2)
p22=p1+p2


time = pseudotime_use+0.1
datt1 = data.frame(time, loc1, loc2)
p1 = ggplot(datt1, aes(x = loc1, y = loc2)) +
       geom_point( alpha =1,size=pointsize,aes(color=clusterlabels)) +
        theme_void()+
        scale_colour_manual(values=color_in)+
    theme(plot.title = element_text(size = textsize,  face = "bold"),
              text = element_text(size = textsize),
              legend.position = "bottom")
datt2 = data.frame(start_x_dat, start_y_dat, end_x_dat, end_y_dat)
p2= geom_segment(aes(x = end_x_dat, y = end_y_dat, xend = start_x_dat, yend = start_y_dat,colour = "segment"),
    arrow = arrow(length = unit(arrowlength,"cm")),size=arrowsize,arrow.fill="black",data = datt2)
p33=p1+p2


return(list("Pseudotime"=p01,"Arrowplot1"=p02,"Arrowplot2"=p03,"Arrowoverlay1"=p22,"Arrowoverlay2"=p33))

}




#' @title Calculate CHAOS score to measure clustering performance.
#' @description CHAOS score measures the spatial continuity of the detected spatial domains.
#' Lower CHAOS score indicates better spatial domian clustering performance.
#' @param clusterlabel Cluster labels.
#' @param location A n by k matrix of spatial locations.
#' @return A numeric value for CHAOS score.
#'
#' @import parallel
#'
#' @export
fx_CHAOS = function(clusterlabel, location){
  # require(parallel)
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  matched_location = scale(matched_location)
  dist_val = rep(0,length(unique(clusterlabel)))
  count = 0
  for(k in unique(clusterlabel)){
    count = count + 1
    location_cluster = matched_location[which(clusterlabel == k),]
    if(length(location_cluster)==2){next}
    #require(parallel)
    results = mclapply(1:dim(location_cluster)[1], fx_1NN, location_in=location_cluster,mc.cores = 5)
    dist_val[count] = sum(unlist(results))
  }
  dist_val = na.omit(dist_val)
  return(sum(dist_val)/length(clusterlabel))

}


#' @title Calculate PAS score to measure clustering performance.
#' @description PAS score measures the randomness of the spots that located outside of the spatial region where it was clustered to.
#' Lower PAS score indicates better spatial domian clustering performance.
#' @param clusterlabel Cluster labels.
#' @param location A n by k matrix of spatial locations.
#' @return A numeric value for PAS score.
#'
#' @import parallel
#'
#' @export
fx_PAS = function(clusterlabel, location){
  # require(parallel)

  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }

  results = mclapply(1:dim(matched_location)[1], fx_kNN, location_in=matched_location,k=10,cluster_in=clusterlabel, mc.cores = 5)
  return(sum(unlist(results))/length(clusterlabel))
}


#' @import pdist
fx_1NN = function(i,location_in){
  # library(pdist)
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  return(min(line_i))
}

#' @import pdist
fx_kNN = function(i,location_in,k,cluster_in){
  #library(pdist)
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  ind = order(line_i)[1:k]
  cluster_use = cluster_in[-i]
  if(sum(cluster_use[ind] != cluster_in[i])>(k/2)){
    return(1)
  }else{
    return(0)
  }

}



# mapDrugToColor<-function(annotations){
#     colorsVector = ifelse(annotations["category"]=="Cluster1",
#         "red", ifelse(annotations["category"]=="Cluster2",
#         "orange",ifelse(annotations["category"]=="Cluster3",
#         "yellow",ifelse(annotations["category"]=="Cluster4",
#         "green",ifelse(annotations["category"]=="Cluster5",
#         "blue",ifelse(annotations["category"]=="Cluster6",
#         "purple",ifelse(annotations["category"]=="Cluster7",
#            "skyblue",ifelse(annotations["category"]=="Cluster8",
#         "black","grey"))))))))
#     return(colorsVector)
# }

# testHeatmap3<-function(logCPM, annotations) {
#     sampleColors = mapDrugToColor(annotations)
#     # Assign just column annotations
#     heatmap3(logCPM, margins=c(10,10),
#       ColSideColors=sampleColors,scale="none",
#       col = colorRampPalette(c( "#0072B2","#F0E442", "#D16103"))(1024),
#       Rowv=NA,
#       Colv=NA,
#         xlab = "Cell ID",
#     ylab = "Marker genes",
#     showColDendro = F,
#   showRowDendro = F)
#     #Assign column annotations and make a custom legend for them
#     heatmap3(logCPM, margins=c(10,10), ColSideColors=sampleColors,
#       scale="none",
#       col = colorRampPalette(c( "#0072B2", "#F0E442","#D16103"))(1024),
#         legendfun=function()showLegend(legend=paste0("Cluster",1:7), col=c("red", "orange", "yellow","green","blue","purple","skyblue"), cex=1),
#               Rowv=NA,
#       Colv=NA,
#         xlab = "Cell ID",
#     ylab = "Marker genes",
#     showColDendro = F,
#   showRowDendro = F)

#     #Assign column annotations as a mini-graph instead of colors,
#     #and use the built-in labeling for them
#     ColSideAnn<-data.frame(Cluster=annotations[["category"]])
#     heatmap3(logCPM, ColSideAnn=ColSideAnn,
#         #ColSideFun=function(x)showAnn(x),
#         margins=c(10,10),
#         ColSideWidth=0.8,
#         Rowv=NA,
#       Colv=NA,
#         xlab = "Cell ID",
#       ylab = "Marker genes",
#       showColDendro = F,
#       showRowDendro = F)
# }




# plot_celltype_barplot_total100 = function(clusternum, celltypes, meta_data_RCTD,method,color_in,textsize=22){

#     if(method == "SpatialPCA"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$SpatialPCA_Louvain==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(meta_data_RCTD)[1]*100,2)
#       }
#     }else if (method == "PCA"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$PCA_Louvain==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(meta_data_RCTD)[1]*100,2)
#       }
#     }else if (method == "NMF"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$NMF_Louvain==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(meta_data_RCTD)[1]*100,2)
#       }
#     }else if (method == "HMRF"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$HMRF==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(meta_data_RCTD)[1]*100,2)
#       }
#     }


#     celltype = names(table(match_type))
#     rownames(percentage) = paste0("Cluster",1:clusternum)
#     colnames(percentage) = names(table(match_type))

#     percentage_vec = c(percentage)
#     cluster = c(rep(c(paste0("Cluster",1:clusternum)),celltypes))
#     celltype = c(rep(celltype,each=clusternum))
#     dat = data.frame(cluster, percentage_vec,celltype)
#     dat$cluster = factor(cluster, level=paste0("Cluster",1:clusternum))

#   p = ggplot(dat, aes(y = percentage_vec,
#              x = factor(cluster), fill = celltype)) +        ## global aes
#   scale_fill_manual(values=color_in)+
#   geom_bar(position="stack", stat="identity",width=0.8,color="grey2") +
#   theme_classic()+xlab("")+ylab("")+
#    theme(plot.title = element_text(size = textsize),
#               text = element_text(size = textsize),
#               #axis.title = element_text(face="bold"),
#               #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
#               #axis.text.x=element_blank(),
#               legend.position = "right")# +


# p

# }


# plot_celltype_barplot_each100 = function(clusternum, celltypes, meta_data_RCTD,method,color_in,textsize=22){

#     if(method == "SpatialPCA"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$SpatialPCA_Louvain==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
#       }
#     }else if (method == "PCA"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$PCA_Louvain==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
#       }
#     }else if (method == "NMF"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$NMF_Louvain==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
#       }
#     }else if (method == "HMRF"){
#       percentage = matrix(0,clusternum,celltypes)
#       for(k in 1:clusternum){
#       metadata_sub = meta_data_RCTD[which(meta_data_RCTD$HMRF==k ),]
#       match_type = metadata_sub$celltype
#       percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
#       }
#     }


#     celltype = names(table(match_type))
#     rownames(percentage) = paste0("Cluster",1:clusternum)
#     percentage_vec = c(percentage)
#     cluster = c(rep(c(paste0("Cluster",1:clusternum)),celltypes))
#     celltype = c(rep(celltype,each=clusternum))
#     dat = data.frame(cluster, percentage_vec,celltype)
#     dat$cluster = factor(cluster, level=paste0("Cluster",1:clusternum))



#   p = ggplot(dat, aes(y = percentage_vec,
#              x = factor(cluster), fill = celltype)) +        ## global aes
#   scale_fill_manual(values=color_in)+
#   geom_bar(position="stack", stat="identity",width=0.8,color="grey2") +
#   theme_classic()+xlab("")+ylab("")+
#    theme(plot.title = element_text(size = textsize),
#               text = element_text(size = textsize),
#               #axis.title = element_text(face="bold"),
#               #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
#               #axis.text.x=element_blank(),
#               legend.position = "right")# +


# p

# }

# #' @title Obtain each cluster label visualized on locations.
# #' @description This function visualizes each cluster label on locations.
# #' @param loc_x: A vector of x coordiantes.
# #' @param loc_y: A vector of y coordinates.
# #' @param clusterlabel_in: A vector of integers, the cluster labels for each cell.
# #' @param pointsize: An integer, the point size of each cell.
# #' @param textsize: An integer, the text size in the legend.
# #' @param title_in: A character string, the title you want to display at the top of the figure.
# #' @param color_in: A vector of colors for each cluster.
# #' @return A list of ggplot objects.
# #' @export
# plot_each_cluster = function(loc_x,loc_y, clusterlabel_in, pointsize=3,text_size=15 ,title_in,color_in){
#   x = loc_x
#   y = loc_y
#   cluster = clusterlabel_in
#   count = 0
#   p=list()
#     for(subcluster in 1:length(unique(cluster))){
#       count = count + 1
#       print(count)
#       clusters = rep("Other_clusters",length(cluster))
#       clusters[which(as.character(cluster)==as.character(subcluster))]=paste0("cluster_",subcluster)
#       datt = data.frame(clusters, x, y)
#       p[[count]] = ggplot(datt, aes(x = x, y = y, color = clusters)) +
#         geom_point( alpha = 0.8,size=pointsize) +
#         scale_colour_manual(values = c("#FFDB6D", "#4E84C4")) +
#           ggtitle(paste0(title_in))+
#         theme_void()+
#         theme(plot.title = element_text(size = text_size),
#               text = element_text(size = text_size),
#               #axis.title = element_text(face="bold"),
#               #axis.text.x=element_text(size = 22) ,
#               legend.position = "bottom")


#   }
# p
# }

# #' @title Obtain louvain clustering cluster labels.
# #' @description This function performs louvain clustering on input low dimensional components.
# #' @param knearest: A vector of integers, number of nearest neighbors for KNN graph construction in louvain clustering.
# #' @param latent_dat: A d by n matrix of low dimensional components.
# #' @return A list of clustering results.
# #' \item{cluster_label}{a list of cluster labels, each corresponds to the vector of nearest neighbors.}
# #' \item{cluster_num}{a list of total number of clusters in each cluster label.}
# #' @export
# louvain_clustering = function(knearest, latent_dat){
#   set.seed(1234)

#   suppressMessages(require(FNN))
#   suppressMessages(require(igraph))

#     PCvalues = latent_dat
#   info.spatial = as.data.frame(t(PCvalues))
#   colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))
#   N_k_nearest=knearest
#   clusterlabel = list()
#   p=list()
#   count = 0
#   max_num = c()
#   for(k_nearest in N_k_nearest){
#       count = count + 1
#       print(count)
#       knn.norm = get.knn(as.matrix(t(PCvalues)), k = k_nearest)
#       knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
#       k=k_nearest), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
#       nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
#       nw.norm = simplify(nw.norm)
#       lc.norm = cluster_louvain(nw.norm)
#       info.spatial$louvain = as.factor(membership(lc.norm))
#       max_num[count] = max(membership(lc.norm))
#       clusterlabel[[count]] = as.character(info.spatial$louvain)
#   }
#   return(list("cluster_label"=clusterlabel, "cluster_num"= max_num))
# }



# #' @title Obtain walktrap clustering cluster labels.
# #' @description This function performs walktrap clustering on input low dimensional components.
# #' @param knearest: A vector of integers, number of nearest neighbors for SNN graph construction in walktrap clustering.
# #' @param latent_dat: A d by n matrix of low dimensional components.
# #' @return A list of clustering results.
# #' \item{cluster_label}{a list of cluster labels, each corresponding to the vector of nearest neighbors}
# #' \item{cluster_num}{a list of total number of clusters in each cluster label.}
# #' @export
# walktrap_clustering = function(knearest, latent_dat){
#   set.seed(1234)

#   suppressMessages(require(bluster))
#   suppressMessages(require(igraph))

#     PCvalues = latent_dat
#   info.spatial = as.data.frame(t(PCvalues))
#   colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))
#   N_k_nearest=knearest
#   clusterlabel = list()
#   p=list()
#   count = 0
#   max_num = c()
#   for(k_nearest in N_k_nearest){
#       count = count + 1
#       print(count)
#       g <- makeSNNGraph(as.matrix(t(PCvalues)),k = k_nearest)
#     clusters <- igraph::cluster_fast_greedy(g)$membership
#     g_walk <- igraph::cluster_walktrap(g)
#     info.spatial$walk = as.factor(membership(g_walk))
#     max_num[count] = max(membership(g_walk))
#       clusterlabel[[count]] = as.character(info.spatial$walk)
#   }
#   return(list("cluster_label"=clusterlabel, "cluster_num"= max_num))
# }






