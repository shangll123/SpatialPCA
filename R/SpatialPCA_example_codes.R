# 
# #############################
# # Pre-Required packages
# #############################
# 
# # for processing the raw data
#   	suppressMessages(require(Seurat))
# 
# # for selecting spatial genes
# 	suppressMessages(require(SPARK))
# 
# # for building kernel
# 	suppressMessages(require(tidyr))
# 	suppressMessages(require(parallel))
# 	suppressMessages(require(MASS))
# 	suppressMessages(require(pdist))
# 	suppressMessages(require(Matrix))
# 	suppressMessages(require(RSpectra))
# 
# # for getting reduced dimension from NMF method
# 	suppressMessages(require(scater))
#   	suppressMessages(require(NMF))
# 
# # for clustering (walktrap method and louvain method)
# 	suppressMessages(require(FNN))
# 	suppressMessages(require(igraph))
# 	suppressMessages(require(bluster))
# 
# # for making figures
#   	suppressMessages(require(ggplot2))
# 
# # for visualize tSNE and UMAP results
# 	suppressMessages(require(Rtsne))
#   	suppressMessages(require(umap))
# 
# # for recording memory usage
# 	suppressMessages(require(peakRAM))
#  
# # for trajectory inference
# 
# 	suppressMessages(require(slingshot))
# 
# # to calculate ARI and NMI
# suppressMessages(require(mclust))
# suppressMessages(require(aricode))
# 
# 
# # functions I haven't updated to the package but used in the analysis:
#  plot_cluster = function(location, clusterlabel, pointsize=3,text_size=15 ,title_in,color_in,legend="none"){
#   cluster = clusterlabel
#   loc_x=location[,1]
#   loc_y=location[,2]
#   datt = data.frame(cluster, loc_x, loc_y)
#   p = ggplot(datt, aes(x = location[,1], y = location[,2], color = cluster)) +
#         geom_point( alpha = 1,size=pointsize) +
#         scale_color_manual(values = color_in)+
#         ggtitle(paste0(title_in))+
#         theme_void()+
#         theme(plot.title = element_text(size = text_size,  face = "bold"),
#               text = element_text(size = text_size),
#               #axis.title = element_text(face="bold"),
#               #axis.text.x=element_text(size = 15) ,
#               legend.position =legend)
#   p
# }
# 
# 
# get_PCA = function(expr,PCnum){
#   
#   n=dim(expr)[2]
#   k = dim(expr)[1]
#   output_sub_mean=matrix(0,k,n)
#   for(i_k in 1:k){
#       output_sub_mean[i_k,]=expr[i_k,]-mean(expr[i_k,])
#   }
#   svd_output_sub_mean=svd(output_sub_mean)
#   A_ini=svd_output_sub_mean$u[,1:PCnum] 
#   Z_pca = t(A_ini) %*% output_sub_mean
#   return(Z_pca)
# }
# 
# get_NMF = function(count, PCnum){
#   suppressMessages(require(scater))
#   suppressMessages(require(NMF))
#   expr = log(count+1)
#   res <- calculateNMF(expr, ncomponents = PCnum)
#   Z_NMF = t(res)
#   return(Z_NMF)
# }
# 
# 
# 
# #############################
# # Data example 1:
# #############################
# 
# 
# i=9 # Here we take the 9th sample as example, in total there are 12 samples (numbered as 1-12)
# clusterNum=c(7,7,7,7,5,5,5,5,7,7,7,7) # each sample has different ground truth cluster number
# 
# suppressMessages(require(SPARK))
# library(Seurat)
# library(peakRAM)
# library(ggplot2)
# library(mclust) # ARI
# library(aricode)# NMI
# 
# # load SpatialPCA updated functions
# path_func="/net/mulan/disk2/shanglu/Projects/SpatialPCA/SpatialPCA_package/SpatialPCA_1.1.0_func/"
# source(paste0(path_func,"SpatialPCA.R"))
# source(paste0(path_func,"SpatialPCA_buildKernel.R"))
# source(paste0(path_func,"SpatialPCA_EstimateLoading.R"))
# source(paste0(path_func,"SpatialPCA_highresolution.R"))
# source(paste0(path_func,"SpatialPCA_SpatialPCs.R"))
# source(paste0(path_func,"SpatialPCA_utilties.R"))
# 
# load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/LIBD/LIBD_sample",i,".RData") )
# xy_coords = as.matrix(xy_coords)
# rownames(xy_coords) = colnames(count_sub)
# LIBD = CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project = "SpatialPCA",gene.type="spatial",sparkversion="spark",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)
# 
# mem <- peakRAM({
#   start_time <- Sys.time()
#   LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
#   LIBD = SpatialPCA_EstimateLoading(LIBD,fast=FALSE,SpatialPCnum=20)
#   LIBD = SpatialPCA_SpatialPCs(LIBD, fast=FALSE)
#   end_time <- Sys.time()
#   T = end_time - start_time
# })
# T
# 
# truth = KRM_manual_layers_sub$layer_guess_reordered[match(colnames(LIBD@normalized_expr),colnames(count_sub))]
# ind_na=which(is.na(truth))
# 
# SpatialPCA_result = list()
# SpatialPCA_result$LIBD  = LIBD
# SpatialPCA_result$SpatialPCs  = LIBD@SpatialPCs
# SpatialPCA_result$normalized_expr  = LIBD@normalized_expr
# SpatialPCA_result$location = LIBD@location
# pred_cluster= walktrap_clustering(clusterNum[i],SpatialPCA_result$SpatialPCs,70 )
# SpatialPCA_result$clusterlabel = pred_cluster
# SpatialPCA_result$clusterlabel_refine=refine_cluster_10x(pred_cluster,SpatialPCA_result$location,shape="hexagon")
# SpatialPCA_result$truth = truth # some spots have NA cluster labels
# SpatialPCA_result$ARI_original = adjustedRandIndex(SpatialPCA_result$clusterlabel[-ind_na],SpatialPCA_result$truth[-ind_na])
# SpatialPCA_result$ARI= adjustedRandIndex(SpatialPCA_result$clusterlabel_refine[-ind_na],SpatialPCA_result$truth[-ind_na])
# SpatialPCA_result$NMI = NMI(as.factor(SpatialPCA_result$clusterlabel_refine[-ind_na]),as.factor(SpatialPCA_result$truth[-ind_na]))
# SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel_refine, SpatialPCA_result$location)
# SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel_refine, SpatialPCA_result$location)
# # save(SpatialPCA_result, file = paste0("LIBD_SpatialPCA_sample",i,"_result.RData"))
# 
# 
# #------------------
# # make figure to check the result
# #------------------
# 
# #------------------
# # first figure for clustering regions
# #------------------
# 
# pdf("Test1.pdf",width=4,height=5)
# cbp=c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91")
# clusterlabel=SpatialPCA_result$clusterlabel_refine
# plot_cluster(legend="none",location=SpatialPCA_result$location,clusterlabel,pointsize=1.5,text_size=20 ,title_in=paste0("SpatialPCA"),color_in=cbp)
# dev.off()
# 
# #------------------
# # second figure for trajectory inference
# #------------------
# 
# library(slingshot)
# sim<- SingleCellExperiment(assays = count_sub)
# reducedDims(sim) <- SimpleList(DRM = t(SpatialPCA_result$SpatialPCs))
# colData(sim)$Walktrap <- factor(SpatialPCA_result$clusterlabel_refine)    
# sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="3" ) 
# # in this data we set white matter region as start cluster, one can change to their preferred start region 
# summary(sim@colData@listData)
# # > summary(sim@colData@listData)
# #                   Length Class              Mode   
# # Walktrap          3639   factor             numeric
# # slingshot         3639   PseudotimeOrdering S4     
# # slingPseudotime_1 3639   -none-             numeric
# 
# sim_SpatialPCA=sim
# # save(sim_SpatialPCA, file = "sim_SpatialPCA.RData")
# pseudotime_traj1 = sim@colData@listData$slingPseudotime_1 # only one trajectory inferred
# clusterlabels = SpatialPCA_result$clusterlabel_refine
# gridnum = 10
# color_in = c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","black")
# pdf("Test2.pdf",width=7,height=5)
# p_traj1 = plot_trajectory(pseudotime_traj1, SpatialPCA_result$location,clusterlabels,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
# print(ggarrange( p_traj1[[4]],p_traj1[[1]],
#           ncol = 2, nrow = 1))
# dev.off()	  
# 
# 
# 
# #------------------
# # make tSNE and UMAP plots
# #------------------
# 
# Z_pca=get_PCA(SpatialPCA_result$normalized_expr,20)
# count_use=count_sub[na.omit(match(rownames(SpatialPCA_result$normalized_expr), rownames(count_sub))),na.omit(match(colnames(SpatialPCA_result$normalized_expr), colnames(count_sub)))]
# Z_NMF=get_NMF(as.matrix(count_use),20)
# 
# library(ggpubr) # for ggarrange function
# 
# set.seed(1234)
# p1 = plot_RGB_tSNE(SpatialPCA_result$location,SpatialPCA_result$SpatialPCs,pointsize=2,textsize=15)
# p2 = plot_RGB_tSNE(SpatialPCA_result$location,Z_pca,pointsize=2,textsize=15)
# p3 = plot_RGB_tSNE(SpatialPCA_result$location,Z_NMF,pointsize=2,textsize=15)
# 
# pdf("Test3.pdf",width=12, height=5)
# ggarrange(p1[[2]], p2[[2]], p3[[2]], 
#           # labels = c("A", "B", "C"),
#           ncol = 3, nrow = 1)
# dev.off()
# 
# 
# set.seed(1234)
# p1 = plot_RGB_UMAP(SpatialPCA_result$location,SpatialPCA_result$SpatialPCs,pointsize=2,textsize=15)
# p2 = plot_RGB_UMAP(SpatialPCA_result$location,Z_pca,pointsize=2,textsize=15)
# p3 = plot_RGB_UMAP(SpatialPCA_result$location,Z_NMF,pointsize=2,textsize=15)
# 
# pdf("Test4.pdf",width=12, height=5)
# ggarrange(p1[[2]], p2[[2]], p3[[2]], 
#           # labels = c("A", "B", "C"),
#           ncol = 3, nrow = 1)
# dev.off()
# 
# 
# #############################
# # Data example 2:
# #############################
# 
# suppressMessages(require(SPARK))
# library(Seurat)
# library(peakRAM)
# library(ggplot2)
# library(mclust) # ARI
# library(aricode)# NMI
# path_func="/net/mulan/disk2/shanglu/Projects/SpatialPCA/SpatialPCA_package/SpatialPCA_1.1.0_func/"
# source(paste0(path_func,"SpatialPCA.R"))
# source(paste0(path_func,"SpatialPCA_buildKernel.R"))
# source(paste0(path_func,"SpatialPCA_EstimateLoading.R"))
# source(paste0(path_func,"SpatialPCA_highresolution.R"))
# source(paste0(path_func,"SpatialPCA_SpatialPCs.R"))
# source(paste0(path_func,"SpatialPCA_utilties.R"))
# 
# 
# # load data
# # remotes::install_github("MarioniLab/DropletUtils") # may need this package to load in data
# library(DropletUtils)
# num=8
# her2stdatanum=34
# load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
# load(paste0("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v2/her2st/data/seu.list.single",num,".rv1.3.RData"))
# H_count = seu.list.single@assays$RNA@counts
# rawcount = H_count[match(rownames(SCTcounts),rownames(H_count)),match(colnames(SCTcounts),colnames(H_count))]
# location=metadata[,5:6]
# location=as.matrix(location)
# 
# # run functions
# ST = CreateSpatialPCAObject(counts=rawcount, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
# mem <- peakRAM({
# start_time <- Sys.time()
# ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
# ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
# ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
# end_time <- Sys.time()
# T = end_time - start_time
# })
# mem
# T
# 
# 
# 
# #------------------
# # clustering
# #------------------
# 
# walktrap_cluster_SpatialPCA = walktrap_clustering(7, ST@SpatialPCs,26)
# cbp_spatialpca <- c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4","chocolate1","lightblue2","#F0E442","red","#CC79A7","mediumpurple","seagreen1")
# pdf("Test5.pdf")
# clusterlabel = walktrap_cluster_SpatialPCA
# loc1 = ST@location[,1]
# loc2 = -ST@location[,2]
# datt = data.frame(clusterlabel, loc1, loc2)
# p = ggplot(datt, aes(x = loc1, y = loc2, color = clusterlabel)) +
#             geom_point( alpha = 1,size=7) +
#             scale_color_manual(values = cbp_spatialpca)+
#             theme_void()+
#             theme(plot.title = element_text(size = 10,  face = "bold"),
#               text = element_text(size = 10),
#               legend.position = "bottom")
# print(p)
# dev.off()
# 
# #------------------
# # high resolution prediction 
# #------------------
# 
# ST = SpatialPCA_highresolution(ST)
# walktrap_SpatialPCA_highresolution = walktrap_clustering(7, ST@highPCs,76)
# cbp_spatialpca <- c(  "plum1", "palegreen4","mediumaquamarine", "chocolate1","#F0E442","dodgerblue","lightblue2","red","#CC79A7","mediumpurple","seagreen1")
# pdf("Test6.pdf",width=5,height=5)
# clusterlabel = walktrap_SpatialPCA_highresolution
# loc1=ST@highPos[,1]
# loc2=-ST@highPos[,2]
# datt = data.frame(clusterlabel, loc1, loc2)
# p = ggplot(datt, aes(x = loc1, y = loc2, color = clusterlabel)) +
#             geom_point( alpha = 1,size=2) +
#             scale_color_manual(values = cbp_spatialpca)+
#             theme_void()+
#             ggtitle(paste0("SpatialPCA high resolution"))+
#             theme(plot.title = element_text(size = 20,  face = "bold"),
#               text = element_text(size = 20),
#               legend.position = "bottom")
# print(p)
# dev.off()
# 
# 
# 
# #############################
# # Data example 3: sample size large, takes ~20G memory and ~30 mins
# #############################
# 
# load(paste0("/net/mulan/disk2/shanglu/Projects/SpatialPCA/data/slideseq/slideseq.rds"))
# slideseq = CreateSpatialPCAObject(counts=sp_count, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5, customGenelist=NULL,min.loctions = 20, min.features=20)
# mem_sparse1 <- peakRAM({
# start_time <- Sys.time()
# slideseq = SpatialPCA_buildKernel(slideseq, kerneltype="gaussian", bandwidthtype="Silverman",bandwidth.set.by.user=NULL,sparseKernel=TRUE,sparseKernel_tol=1e-20,sparseKernel_ncore=10)
# slideseq = SpatialPCA_EstimateLoading(slideseq,fast=TRUE,SpatialPCnum=20)
# slideseq = SpatialPCA_SpatialPCs(slideseq, fast=TRUE)
# end_time <- Sys.time()
# T_sparse1 = end_time - start_time
# })
# mem_sparse1
# T_sparse1
# 
# 
# louvain_cluster_SpatialPCA = louvain_clustering( 8, as.matrix(slideseq@SpatialPCs),260)
# clusterlabel = as.character(louvain_cluster_SpatialPCA)
# cbp_spatialpca <- c("#66C2A5", "lightyellow2", "cornflowerblue" ,"#E78AC3", "skyblue1" ,"#FFD92F" ,"lightcyan2", "coral")
# 
# pdf("Test7.pdf",width=5,height=5)
# loc1 = slideseq@location[,1]
# loc2 = slideseq@location[,2]
# datt = data.frame(clusterlabel, loc1, loc2)
# p = ggplot(datt, aes(x = loc1, y = loc2, color = clusterlabel)) +
#             geom_point( alpha = 1,size=0.8) +
#             scale_color_manual(values = cbp_spatialpca)+
#             theme_void()+
#             theme(plot.title = element_text(size = 10,  face = "bold"),
#               text = element_text(size = 10),
#               legend.position = "bottom")
# print(p)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
