### dir_current <- "F:/Research paper/IntegrateDRcluster/AnalysisCode/PRECAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/utility_funcs.R"))

setwd(paste0(dir_current, "Real_data_results/dataFiles/HCC4/") )


##### 5b #####

library(cowplot)
load("HCC4_posList.rds")
load("HCC4_RGB_hZ_umapList_allsample.rds")
posList1 <- posList

load("RGB_umap_tsne_paste_hcc4.rds")
str(hZ_umap3List_allsample)
## Add the results of PASTE
hZ_umap3List_allsample[["PASTE"]] <- hZ_umap3_paste

### RGB plot for each sample based Joint UMAP
indexList <- get_indexList(posList)
hZ_umap3List_allsample[2:3] <- hZ_umap3List_allsample[3:2]
names(hZ_umap3List_allsample)[2:3] <- names(hZ_umap3List_allsample)[3:2]


### Plot iDR-SC for 4 samples

pList_idrsc_umap_RGB <- list()
datList <- list()
for(r in 1:4){## each sample
  # i <- 1
  
  pList_eachmethod <- list()
  for(i in 1){  ## each method
    # r <- 1
    message("r = ", r)
    pt_size = 0.5
    if(r == 1) pt_size= 0.3
    dat_tmp <- as.data.frame(cbind(posList1[[r]], hZ_umap3List_allsample[[i]][indexList[[r]],]))
    dat_tmp$Sample <- paste0("HCC", r)
    colnames(dat_tmp) <- c("Coord x", "Coord y", "UMAP1", "UMAP2", "UMAP3", "Sample")
    datList[[r]] <- dat_tmp
    ptmp <- plot_RGB(posList1[[r]], hZ_umap3List_allsample[[i]][indexList[[r]],], pointsize = 0.8) +
      mytheme_graybox(border_color = 'white') + xlim(c(0,77)) + ylim(c(1,120)) + geom_point(alpha=0.5)
    pList_idrsc_umap_RGB[[r]]<- ptmp 
  }
  
}

p12 <- plot_grid(plotlist =pList_idrsc_umap_RGB, nrow=1, ncol=4)
ggsave(file="HCC4_umap_RGB_idrsc.png", plot = p12,
       width = 14, height =4,  units = "in", dpi = 500)



load("cluster_SCMEB_methods7_HCC4.rds")
load('HCC4_clusterK9_renumber_idrsc.rds')
load("cluster_allsample_paste_hcc4.rds")
str(yMat_paste)
clusterMat <- cbind("iDR-SC"=unlist(cluster_idrsc_renumber), clusterMat, "PASTE" = yMat_paste[,1])
head(clusterMat)

MethodNames <- c("iDR-SC", "Seurat V3" ,"Harm-SC-MEB",
                 "fastMNN-SC-MEB",  "Scanorama-SC-MEB",
                 "scGen-SC-MEB", "scVI-SC-MEB", 'MEFISTO-SC-MEB', "PASTE")
clusterMat_sub <- clusterMat[, MethodNames]
dim(clusterMat_sub)
apply(clusterMat_sub, 2, max)
clusterMat_sub <- matrix(as.numeric(clusterMat_sub), ncol= length(MethodNames))

load("HCC4_posList.rds")
apply(clusterMat, 2, max)
indexList <- get_indexList(posList)
cols_cluster <- c("#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1",
                  "#1f77b4", "#A65628",  "blue", "green")
### plot all samples for PRECAST
pList_assign_iDRSC <- list()
for(j in 1){ ## each method
  # j <- 1
  pList_eachmethod <- list()
  for(r in 1:4){ ## each sample
    message("r = ", r)
    cluster_tmp <- factor(clusterMat_sub[indexList[[r]],j])
    datList[[r]]$Domain <- cluster_tmp
    p_tmp <- plot_scatter(posList1[[r]], meta_data = data.frame(cluster=cluster_tmp),
                          label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                          point_size = 1.0,palette_use = cols_cluster[as.numeric(levels(cluster_tmp))], point_alpha = 0.6) + 
      mytheme_graybox(border_color = 'white') + xlim(c(0,77)) + ylim(c(1,120))
    
    pList_assign_iDRSC[[r]]<- p_tmp
  }
}

p12 <-cowplot:: plot_grid(plotlist =pList_assign_iDRSC, nrow=1, ncol=4)
ggsave(file="HCC4_assign_heatmap_iDRSC.png", plot = p12,
       width = 14, height =4, units = "in", dpi = 500)



##### 5c #####

##tSNE plot VS sample and clusters
load("hZtsneList_HCC4.rds")
load('HCC4_clusterK9_renumber_idrsc.rds')
load("tsne2_paste_hcc4.rds")
hZtsneList[["PASTE"]] <- hZ_tsne2_paste

head(meta_data)
meta_data$cluster_reorder <- factor(unlist(cluster_idrsc_renumber))
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names(palettes)
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_pal <- c(pal1(max_n), pal2(max_n2)[c(1,3,8,12)]);
# cols_cluster <- cols_pal
# simutool::colorbar_adj_transparent(cols_cluster[1:9])
library(RColorBrewer)
RColorBrewer::display.brewer.all()

cols_cluster <- c("#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1",
                  "#1f77b4", "#A65628",  "blue", "green")

require(colorspace)
ramp.list = adjust_transparency("green",   alpha = 0.1)
cols_cluster2 <- cols_cluster
cols_cluster2[9] <- 'green3'
## select colors

cols_sample <- gg_color_hue(4)


pt_size_sample <- 0.8; pt_alpha_sample <- 0.5
pt_size_cluster <- 0.8; pt_alpha_cluster <- 1
base_axis_size <- 20

library(ggplot2)
str(hZtsneList)
pList_sample <- list(); pList_cluster <- list()
for(j in 1:length(hZtsneList)){
  message('j = ', j)
  
  
  pList_sample[[j]] <- plot_scatter(hZtsneList[[j]], meta_data, label_name = 'sample',border_col="gray10", 
                                    base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
                                    point_alpha = pt_alpha_sample,no_guides=T) + labs(x=NULL, y=NULL)
  pList_cluster[[j]] <- plot_scatter(hZtsneList[[j]], meta_data, label_name = 'cluster_reorder',
                                     base_size=20,border_col="gray10",
                                     palette_use=cols_cluster2,
                                     point_size = pt_size_cluster, 
                                     point_alpha =pt_alpha_cluster,no_guides=T) + labs(x=NULL, y=NULL)
  if(j==2 || j==5){
    pList_sample[[j]] <- pList_sample[[j]] +xlim(c(-40, 40))
    pList_cluster[[j]] <- pList_cluster[[j]] +xlim(c(-40, 40))
  }
  
}
names(pList_sample) <-names(pList_cluster) <-  names(hZtsneList)

library(cowplot)
#pList_sample[[1]]
subnames <- c("iDR-SC" ,"Seurat V3", "Harmony","fastMNN", "Uncorrected")
p12 <- plot_grid(plotlist =pList_sample[subnames], nrow=1, ncol=5)

ggsave(file="HCC4_tsne_sample_heatmap.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 500)
p1 <- pList_cluster[[1]]
pList_cluster[[1]] <- Seurat::  LabelClusters(p1, id='cluster_reorder', repel = T, 
                                              size = 9, color = 'white', 
                                              box = T, position = "nearest")
p12 <- plot_grid(plotlist =pList_cluster[subnames], nrow=1, ncol=5)
ggsave(file="HCC4_tsne_cluster_heatmap.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 500)


# subnames2 <- c("Scanorama",  "scGen", "scVI", "MEFISTO", "PASTE")
# p12 <- plot_grid(plotlist =pList_sample[subnames2], nrow=1, ncol=5)
# ggsave(file="HCC4_tsne_sample_heatmap_Py4V2.png", plot = p12,
#        width = 25, height =5, units = "in", dpi = 200)
# p12 <- plot_grid(plotlist =pList_cluster[subnames2], nrow=1, ncol=5)
# ggsave(file="HCC4_tsne_cluster_heatmap_Py4V2.png", plot = p12,
#        width = 25, height =5, units = "in", dpi = 200)


plot_scatter(hZtsneList[[1]], meta_data, label_name = 'cluster_reorder', base_size=20,
             palette_use=cols_cluster2, 
             point_size = pt_size_cluster, 
             point_alpha =pt_alpha_cluster,no_guides=F)


##### 5d & 5e #####
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggthemes)
library(cowplot)
library(Matrix)
library(data.table)
library(reshape2)

load("HCC4_clusterK9_renumber_idrsc.rds")
load("HCC4_posList.rds")


color_pal = list("#CB181D","#EF3B2C","#FB6A4A","#FC9272","#FCBBA1","#1f77b4","#ff7f0e","#2ca02c","#c5b0d5")
slice_set_sequence = seq(1,4)
need_scale = T
datList <- list()
figure_list = list()
for (slice_set in slice_set_sequence){
  norm_weights = read.csv(paste('hcc_deconvolution/output_weights_geo125set1_rctd_st',slice_set,'.csv',sep=""),row.names=1)
  norm_weights$Immune = norm_weights$TEC+norm_weights$T.cell+norm_weights$B.cell
  norm_weights = norm_weights[,c('Immune',"CAF","TAM","HPC.like",'Malignant.cell')]
  
  cluster_id = cbind(posList[[slice_set]],cluster_idrsc_renumber[[slice_set]])
  cluster_weight = merge(norm_weights,cluster_id[,3],by = 0)
  rownames(cluster_weight) = cluster_weight$Row.names
  cluster_weight = cluster_weight[,-1]
  
  percentage = as.data.frame(cluster_weight %>% group_by(y) %>% summarise(across(everything(), sum)))
  colnames(percentage) = c( "Cluster",'Immune cell',"CAF","TAM","HPC-like cell","Malignant cell")
  percentage_long <- melt(as.data.frame(percentage),id.vars ='Cluster')
  
  if (need_scale == T){
    percentage_scale = matrix(0,nrow(percentage),ncol(percentage))
    for (tt in 1:nrow(percentage)){
      percentage_scale[tt,] = as.numeric(cbind(percentage$Cluster[tt],percentage[tt,2:ncol(percentage)]/sum(rowSums(percentage[,2:ncol(percentage)]))))
    }
    colnames(percentage_scale) = colnames(percentage)
    #sum(percentage_scale[,2:ncol(percentage)])
    percentage_long <- melt(as.data.frame(percentage_scale),id.vars ='Cluster')
  }
  
  if (need_scale == T){
    geom_bar_position = 'stack'
  }else{
    geom_bar_position = 'fill'
  }
  level_names <- c("Malignant cell","Immune cell","CAF","TAM","HPC-like cell") 
  percentage_long$variable <- factor(percentage_long$variable, levels=level_names)
  levels(percentage_long$variable)
  percentage_long$Sample <- paste0("HCC", slice_set)
  datList[[slice_set]] <- percentage_long
  if (slice_set == 4){
    figure_list[[slice_set]] <- bar_plot(percentage_long,geom_bar_position=geom_bar_position, legend_position='right',
                                         color_pal = ggthemes::tableau_color_pal()(5)[c(3,2,5,1,4)])}else{
                                           figure_list[[slice_set]] <- bar_plot(percentage_long,geom_bar_position=geom_bar_position, legend_position='none',
                                                                                color_pal = ggthemes::tableau_color_pal()(5)[c(3, 2,5,1,4)])
                                         }
}

figure_celltype_weight <- ggarrange(plotlist = figure_list,ncol=2,nrow=2,align='hv',common.legend = T,
                                    legend = 'right')
if(need_scale){
  ggsave(file='result/hcc_percentage_scaled.png', plot = figure_celltype_weight, 
         width = 10, height = 8, units = "in", bg = 'white', dpi = 500,limitsize = F)
}else{
  ggsave(file='result/hcc_percentage_unscaled.png', plot = figure_celltype_weight, 
         width = 10, height = 8, units = "in", bg = 'white', dpi = 500,limitsize = F)
}



load("HCC4_clusterK9_renumber_idrsc.rds")
load("HCC4_posList.rds")

color_pal = list("#CB181D","#EF3B2C","#FB6A4A","#FC9272","#FCBBA1","#1f77b4","#ff7f0e","#2ca02c","#c5b0d5")

cols_cluster <- c("#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1",
                  "#1f77b4", "#A65628",  "blue", "green2")
# install.packages("pals")
# immune
my_pal = pals::brewer.greys(20)[2:20]
figure_list_immune = list()
for (slice_set in 1:4){
  # slice_set <- 1
  
  st_coord <- as.data.frame(posList[[slice_set]])
  colnames(st_coord) <- c("imagerow", "imagecol")
  
  norm_weights = read.csv(paste('hcc_deconvolution/output_weights_geo125set1_rctd_st',slice_set,'.csv',sep=""),row.names=1)
  norm_weights$Immune = norm_weights$TEC+norm_weights$T.cell+norm_weights$B.cell
  norm_weights = norm_weights[,c("CAF","TAM","HPC.like",'Malignant.cell','Immune')]
  
  cluster_id = cbind(posList[[slice_set]],cluster_idrsc_renumber[[slice_set]])
  cluster_weight = merge(norm_weights,cluster_id[,3],by = 0)
  rownames(cluster_weight) = cluster_weight$Row.names
  cluster_weight = cluster_weight[,-1]
  
  #st_coord$imagerow = max(st_coord$imagerow)-st_coord$imagerow
  
  ct=5
  plot_val = as.data.frame(cluster_weight[,ct])
  rownames(plot_val) = rownames(cluster_weight)
  colnames(plot_val) = colnames(cluster_weight)[ct]
  barcodes = rownames(plot_val)
  my_table = st_coord[barcodes, ]
  my_table$Proportion = plot_val[barcodes,]
  my_table$Proportion <- range01(my_table$Proportion)
  med <- quantile(my_table$Proportion, 0.94)
  my_table$Proportion[ my_table$Proportion > med] <- med
  
  my_table$cluster = cluster_weight[barcodes,'y']
  ylimit = c(0, 1)
  
  my_table$border_color = '#cccccc'
  my_table$stroke = 0
  
  select_cluster = c(6,7,8,9)
  for (i in select_cluster){
    my_table[my_table$cluster==i,'stroke'] = 0.45
  }
  
  for (i in 1:10){
    my_table[my_table$cluster==i,"border_color"] = cols_cluster[i]#color_pal[i]
  }
  
  table(my_table$cluster)
  
  plot <- ggplot(my_table, aes(x =imagerow , y = imagecol)) + 
    geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 5,pch=21) +
    #scale_fill_gradientn(colors = my_pal,limits = ylimit) + ,colour = my_table$border_color
    scale_fill_gradientn(colors = c("#361A95","white", "#D62728")) +
    scale_shape_identity() + scale_size_identity() + theme_classic() +
    mytheme_graybox(base_size = 28)
  
  
  figure_list_immune[[slice_set]] <- plot +xlim(c(0,77)) + ylim(c(1,120))# + coord_fixed() 
  if (slice_set == 1){
    figure_list_immune[[slice_set]] <- figure_list_immune[[slice_set]] + ggtitle("Immune cell")
  }
}
figure_celltype_weight <- ggarrange(plotlist = figure_list_immune,ncol=4,nrow=1,align='hv',common.legend = T,
                                    legend = 'right')
figure_celltype_weight
ggsave(file='result/Immune_select_type_renumber.png', plot = figure_celltype_weight, 
       width = 24, height = 6, units = "in", bg = 'white', dpi = 500,limitsize = F)



# Malig
my_pal = pals::brewer.greys(20)[2:20]
figure_list_malig = list()
for (slice_set in 1:4){
  st_coord <- as.data.frame(posList[[slice_set]])
  colnames(st_coord) <- c("imagerow", "imagecol")
  
  
  norm_weights = read.csv(paste('hcc_deconvolution/output_weights_geo125set1_rctd_st',slice_set,'.csv',sep=""),row.names=1)
  norm_weights$Immune = norm_weights$TEC+norm_weights$T.cell+norm_weights$B.cell
  norm_weights = norm_weights[,c("CAF","TAM","HPC.like",'Malignant.cell','Immune')]
  
  cluster_id = cbind(posList[[slice_set]],cluster_idrsc_renumber[[slice_set]])
  cluster_weight = merge(norm_weights,cluster_id[,3],by = 0)
  rownames(cluster_weight) = cluster_weight$Row.names
  cluster_weight = cluster_weight[,-1]
  
  
  
  ct=4
  plot_val = as.data.frame(cluster_weight[,ct])
  rownames(plot_val) = rownames(cluster_weight)
  colnames(plot_val) = colnames(cluster_weight)[ct]
  barcodes = rownames(plot_val)
  my_table = st_coord[barcodes, ]
  my_table$Proportion = plot_val[barcodes,]
  
  my_table$Proportion <- range01(my_table$Proportion)
  med <- quantile(my_table$Proportion, 0.94)
  my_table$Proportion[ my_table$Proportion > med] <- med
  
  
  my_table$cluster = cluster_weight[barcodes,'y']
  ylimit = c(0, 1)
  
  my_table$border_color = '#cccccc'
  my_table$stroke = 0
  
  select_cluster = c(6,7,8,9)
  for (i in select_cluster){
    my_table[my_table$cluster==i,'stroke'] = 0.45
  }
  for (i in 1:10){
    my_table[my_table$cluster==i,"border_color"] = cols_cluster[i]#ggthemes::tableau_color_pal()(10)[i]
  }
  
  table(my_table$cluster)
  head(my_table)
  
  plot <- ggplot(my_table, aes(x = imagerow, y = imagecol)) + 
    geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 3,pch=21) +
    #scale_fill_gradientn(colors = my_pal,limits = ylimit) + colour = my_table$border_color
    scale_fill_gradientn(colors = c("#361A95","white", "#D62728")) +
    scale_shape_identity() + scale_size_identity() + theme_classic() +
    mytheme_graybox(base_size = 28)
  
  # xlim <- c(min(st_coord$imagecol) - 1, max(st_coord$imagecol) + 1)
  # ylim <- c(min(st_coord$imagerow) - 1, max(st_coord$imagerow) + 1)
  figure_list_malig[[slice_set]] <- plot  +xlim(c(0,77)) + ylim(c(1,120))#+ coord_fixed() + xlim(xlim) + ylim(ylim)
  if (slice_set == 1){
    figure_list_malig[[slice_set]] <- figure_list_malig[[slice_set]] + ggtitle("Malignant cell")
  }
}
figure_celltype_weight <- ggarrange(plotlist = figure_list_malig,ncol=4,nrow=1,
                                    align='hv',common.legend = T,
                                    legend = 'right')

ggsave(file='result/Malig_select_type.png', plot = figure_celltype_weight, 
       width = 24, height = 6, units = "in", bg = 'white', dpi = 500,limitsize = F)



# HPC.like
my_pal = pals::brewer.greys(20)[2:20]
figure_list_hpc = list()
for (slice_set in 1:4){
  
  
  st_coord <- as.data.frame(posList[[slice_set]])
  colnames(st_coord) <- c("imagerow", "imagecol")
  
  norm_weights = read.csv(paste('hcc_deconvolution/output_weights_geo125set1_rctd_st',slice_set,'.csv',sep=""),row.names=1)
  norm_weights$Immune = norm_weights$TEC+norm_weights$T.cell+norm_weights$B.cell
  norm_weights = norm_weights[,c("CAF","TAM","HPC.like",'Malignant.cell','Immune')]
  
  cluster_id = cbind(posList[[slice_set]],cluster_idrsc_renumber[[slice_set]])
  cluster_weight = merge(norm_weights,cluster_id[,3],by = 0)
  rownames(cluster_weight) = cluster_weight$Row.names
  cluster_weight = cluster_weight[,-1]
  
  # st_coord$imagerow = max(st_coord$imagerow)-st_coord$imagerow
  
  ct=3
  plot_val = as.data.frame(cluster_weight[,ct])
  rownames(plot_val) = rownames(cluster_weight)
  colnames(plot_val) = colnames(cluster_weight)[ct]
  barcodes = rownames(plot_val)
  my_table = st_coord[barcodes, ]
  my_table$Proportion = plot_val[barcodes,]
  my_table$Proportion <- range01(my_table$Proportion)
  med <- quantile(my_table$Proportion, 0.94)
  my_table$Proportion[ my_table$Proportion > med] <- med
  
  my_table$cluster = cluster_weight[barcodes,'y']
  ylimit = c(0, 1)
  
  my_table$border_color = '#cccccc'
  my_table$stroke = 0
  
  select_cluster = c(6,7,8,9)
  for (i in select_cluster){
    my_table[my_table$cluster==i,'stroke'] = 0.42
  }
  
  for (i in 1:10){
    my_table[my_table$cluster==i,"border_color"] = cols_cluster[i] #ggthemes::tableau_color_pal()(10)[i]
  }
  
  table(my_table$cluster)
  head(my_table)
  
  plot <- ggplot(my_table, aes(x = imagerow, y = imagecol)) + 
    geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 3.5,pch=21) +
    #scale_fill_gradientn(colors = my_pal,limits = ylimit) +  # ,colour = my_table$border_color
    scale_fill_gradientn(colors = c("#361A95","white", "#D62728")) +
    scale_shape_identity() + scale_size_identity() + theme_classic() +
    mytheme_graybox(base_size = 28)
  
  xlim <- c(min(st_coord$imagecol) - 1, max(st_coord$imagecol) + 1)
  ylim <- c(min(st_coord$imagerow) - 1, max(st_coord$imagerow) + 1)
  figure_list_hpc[[slice_set]] <- plot  +xlim(c(0,77)) + ylim(c(1,120)) #+ coord_fixed() + xlim(xlim) + ylim(ylim)
  if (slice_set == 1){
    figure_list_hpc[[slice_set]] <- figure_list_hpc[[slice_set]] + ggtitle("HPC-like cell")
  }
}
figure_celltype_weight <- ggarrange(plotlist = figure_list_hpc,ncol=4,nrow=1,align='hv',common.legend = T,
                                    legend = 'right')

ggsave(file='result/HPC_select_type7.png', plot = figure_celltype_weight, 
       width = 24, height = 6, units = "in", bg = 'white', dpi = 500,limitsize = F)




