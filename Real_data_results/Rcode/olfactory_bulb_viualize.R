### dir_current <- "F:/Research paper/IntegrateDRcluster/AnalysisCode/PRECAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/utility_funcs.R"))
setwd(paste0(dir_current, "Real_data_results/dataFiles/Bulb16/") )


##### 4b #####
#### Start Assignment plot
load("clusterAll_tmpList_Method8_Bulb16_rep2_merge70.rds")
load('posList16_before_merge70_Bulb16rep2.rds')
posList <- posList16_before_merge70
sum(sapply(posList, nrow))
## Add PASTE
load("clusterAll_paste_tmpList_Bulb16_rep2_merge70.rds")
for(r in 1:16){
  clusterAll_tmpList[[r]][["PASTE"]] <- clusterAll_paste_tmpList[[r]]
}

idx_outList <- list()
for(r in 3:16){
  idx_outList[[r]] <-  c(order(posList[[r]][,2], decreasing = T)[1:10],
                         order(posList[[r]][,1])[1:10])
  
}
##### domain color selection
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names(palettes)
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")


cols_cluster <- c( "#FD7446", 'blue3',"#709AE1", "green1",
                   "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
                   "#91D1C2","#FFD70099", "#6B6ECF",  "#C7E9C0" ,
                   "#7B4173", "green4","#FF9896",  pal1(max_n))

## change the order of clusters
rawID <-  c(1, 5, 3, 12, 9,11, 7,  8, 4, 10, 2,  6)
names(rawID) <- 1: 12
# simutool::colorbar_adj_transparent(cols_cluster[1:12], alpha=1)
#####plot iDR-SC cluster

pList <- list()
pt_size <- 0.2
point_alpha <- 0.8
datList <- list()
for(r in 1:16){
  # r <- 8
  message('r = ', r)
  if(r >=3){
    cluster_tmp <- clusterAll_tmpList[[r]][[1]][-idx_outList[[r]]]
  }else{
    cluster_tmp <- clusterAll_tmpList[[r]][[1]]
  }
  
  idx <- !is.na(cluster_tmp)
  cluster_tmp <- as.numeric(replace_ID(cluster_tmp[idx], rawID))
  unique_sort <- sort(unique(cluster_tmp))
  clusterlabel <- factor(cluster_tmp, levels=unique_sort)
  pos_tmp <- if(r>=3){
    posList[[r]][-idx_outList[[r]],]
  }else{
    posList[[r]] 
  } 
  if(r != 9){
    datList[[r]] <- cbind(pos_tmp[idx, ], clusterlabel)
    ptmp <- plot_scatter(pos_tmp[idx, ], 
                         meta_data = data.frame(cluster=clusterlabel),
                         label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                         point_size = pt_size,palette_use = cols_cluster[unique_sort], point_alpha = point_alpha) + 
      theme_bulb16()
  }else{
    datList[[r]] <- cbind(pos_tmp[idx, ][-32115 ,], clusterlabel)
    ptmp <- plot_scatter(pos_tmp[idx, ][-32115 ,], meta_data = data.frame(cluster=clusterlabel),
                         label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                         point_size = pt_size,palette_use = cols_cluster[unique_sort], point_alpha = point_alpha) + 
      theme_bulb16()
  }
  
  
  pList[[r]] <- ptmp
}


p12 <- cowplot::plot_grid(plotlist=pList, nrow=3, ncol=6)
# ggsave(file="Bulb16_rep2_merge70_mPotts5_chooseK12_clusterAssignment.png", plot = p12,
#        width = 11, height =11, units = "in", dpi = 500)

ggsave(file="Bulb16_rep2_merge70_mPotts5_chooseK12_clusterAssignmentV2.png", plot = p12,
       width = 12, height =6, units = "in", dpi = 500)

## Label legend
p1 <- pList[[8]]
Seurat:: LabelClusters(p1, id='cluster', repel = T, 
                       size = 5, color = 'white', 
                       box = T, position = "nearest", max.overlaps =8)

plot_scatter(pos_tmp[idx, ], meta_data = data.frame(cluster=clusterlabel),
             label_name = 'cluster', xy_names = c("", ""), no_guides = F,
             point_size = pt_size,palette_use = cols_cluster[unique_sort], point_alpha = point_alpha) + 
  mytheme_graybox()

##### 4c #####
### tSNE plots vs sample and clusters
load("tsne2_allMethods_Bulb16_rep2_merge70.rds")

load("tsne2_paste_bulb16_merge.rds")
tsne2List[["PASTE"]] <- hZ_tsne2_paste

meta_data$domain <- factor(replace_ID(meta_data$domain, rawID), levels=1:12)
cols_sample=c('#1CE6FF','#FF34FF','#008941','#A30059',
              '#FF2F80','#0000A6','#63FFAC','#004D43','#8FB0FF','#4FC601',
              '#3B5DFF','#4A3B53','#61615A','#BA0900','#6B7900','#00C2A0',
              '#FFAA92','#FF90C9','#B903AA')
# colorbar_adj_transparent(cols_cluster[1:12],1)

cols_cluster2 <- cols_cluster
cols_cluster2[4] <- "green2"
library(ggplot2)
str(tsne2List)


pt_size_sample <- 0.2; pt_alpha_sample <- 0.5
pt_size_cluster <- 0.2; pt_alpha_cluster <- 0.7
base_axis_size <- 20

pList_sample <- list(); pList_cluster <- list()
for(j in 1:length(tsne2List)){
  # j <- 1
  message('j = ', j)
  if(j %in% c(1,4:10)){
    pList_sample[[j]] <- plot_scatter(tsne2List[[j]], meta_data, label_name = 'sample', border_col='gray13', 
                                      base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
                                      point_alpha = pt_alpha_sample,no_guides=T) + labs(x=NULL, y=NULL)
    pList_cluster[[j]] <- plot_scatter(tsne2List[[j]], meta_data, label_name = 'domain', base_size=20,
                                       palette_use=cols_cluster2, point_size = pt_size_cluster,border_col='gray13', 
                                       point_alpha =pt_alpha_cluster,no_guides=T) + labs(x=NULL, y=NULL)
    
    if(j %in% c(2,3)){
      pList_sample[[j]] <- pList_sample[[j]] + xlim(c(-40, 30)) 
      pList_cluster[[j]] <- pList_cluster[[j]] + xlim(c(-40, 30)) 
    }
  }
}
names(pList_sample) <-names(pList_cluster) <-  names(tsne2List)

library(cowplot)
subnames <- c("iDR-SC" ,"Seurat V3", "Harmony","fastMNN", "Uncorrected")
p12 <- plot_grid(plotlist =pList_sample[subnames], nrow=1, ncol=5)
ggsave(file="tsne_sample_heatmap_Bulb16_rep2_merge70.png", plot = p12,
       width = 25, height =4.8, units = "in", dpi = 500)

p1 <- pList_cluster[[1]]
pList_cluster[[1]] <- Seurat::  LabelClusters(p1, id='domain', repel = T, 
                                              size = 9, color = 'white', 
                                              box = T, position = "nearest")

p12 <- plot_grid(plotlist =pList_cluster[subnames], nrow=1, ncol=5)
ggsave(file="tsne_cluster_heatmap_Bulb16_rep2_merge70.png", plot = p12,
       width = 25, height =4.8, units = "in", dpi = 500)

## get legend Label
j <- 1
plot_scatter(tsne2List[[j]], meta_data, label_name = 'sample', border_col='gray13', 
             base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
             point_alpha = pt_alpha_sample,no_guides=F)

plot_scatter(tsne2List[[j]], meta_data, label_name = 'domain', base_size=20,
             palette_use=cols_cluster2, point_size = pt_size_cluster,border_col='gray13', 
             point_alpha =pt_alpha_cluster,no_guides=F)

# pList_sample[[5]]
# subnames2 <- c("Scanorama",  "scGen", "scVI", "MEFISTO", "PASTE")
# p12 <- plot_grid(plotlist =pList_sample[subnames2], nrow=1, ncol=5)
# ggsave(file="tsne_Py4_sample_heatmap_Bulb16_rep2_merge70.png", plot = p12,
#        width = 25, height =5, units = "in", dpi = 200)
# p12 <- plot_grid(plotlist =pList_cluster[subnames2], nrow=1, ncol=5)
# ggsave(file="tsne_Py4_cluster_heatmap_Bulb16_rep2_merge70.png", plot = p12,
#        width = 25, height =5, units = "in", dpi = 200)

##### 4d #####
library(ggthemes)
cols_cluster <- c(get_trans_colors("gray16", 5)[1:3], "#1f77b4", get_trans_colors("#bcbd22", 4)[2:3],
                  "darkolivegreen4",
                  get_trans_colors("#709AE1",7),"#ff7f0e", get_trans_colors("#ffbb78", 3), "darkorchid2",
                  get_trans_colors("#31A354",2), 
                  get_trans_colors("#DE9ED6", 3),"#BCBD22", c('deepskyblue2',
                                                              'deepskyblue3'),#get_trans_colors("deepskyblue2",2),
                  "#DADAEB",
                  get_trans_colors("#91D1C2",5),  "#6B6ECF", "#7B4173", 
                  get_trans_colors("#9EDAE5", 3), "#d62728", "#ff9896", "red", 'red1' )


#########bar plot for

load("Bulb16_rep2_Merge70_clusterK12_renumber_idrsc.rds")
table(cluster_idrsc_renumber[[1]])
load("posList_Bulb16_rep2_merge70.rds")
str(posList[[1]])


slice_set_sequence = 1
need_scale = T
figure_list = list()
datList <- list()
for (slice_set in slice_set_sequence){
  # slice_set <- 1
  norm_weights = read.csv(paste('deconvoltion_results_new/output_weights_bulb_st',slice_set,'.csv',sep=""),row.names=1)
  
  colnames(norm_weights)[15] <- "Macrophage" #'Mf'
  
  
  colnames(norm_weights)[colnames(norm_weights) =='n01.OSNs'] = 'OSNs'
  colnames(norm_weights)[colnames(norm_weights) =='n02.PGC.1'] = 'PGC-1'
  colnames(norm_weights)[colnames(norm_weights) =='n03.GC.1'] = 'GC-1'
  colnames(norm_weights)[colnames(norm_weights) =='n04.Immature'] = 'Immature'
  colnames(norm_weights)[colnames(norm_weights) =='n05.PGC.2'] = 'PGC-2'
  colnames(norm_weights)[colnames(norm_weights) =='n06.Transition'] = 'Transition'
  colnames(norm_weights)[colnames(norm_weights) =='n07.GC.2'] = 'GC-2'
  colnames(norm_weights)[colnames(norm_weights) =='n08.PGC.3'] = 'PGC-3'
  colnames(norm_weights)[colnames(norm_weights) =='n09.GC.3'] = 'GC-3'
  colnames(norm_weights)[colnames(norm_weights) =='n10.GC.4'] = 'GC-4'
  colnames(norm_weights)[colnames(norm_weights) =='n11.GC.5'] = 'GC-5'
  colnames(norm_weights)[colnames(norm_weights) =='n12.GC.6'] = 'GC-6'
  colnames(norm_weights)[colnames(norm_weights) =='n13.AstrocyteLike'] = 'AstrocyteLike'
  colnames(norm_weights)[colnames(norm_weights) =='n14.GC.7'] = 'GC-7'
  colnames(norm_weights)[colnames(norm_weights) =='n15.M.TC.1'] = 'M/TC-1'
  colnames(norm_weights)[colnames(norm_weights) =='n16.M.TC.2'] = 'M/TC-2'
  colnames(norm_weights)[colnames(norm_weights) =='n17.M.TC.3'] = 'M/TC-3'
  colnames(norm_weights)[colnames(norm_weights) =='n18.EPL.IN'] = 'EPL-IN'
  
  
  cluster_id = cbind(posList[[slice_set]],cluster_idrsc_renumber[[slice_set]])
  cluster_weight = merge(norm_weights,cluster_id[,3],by = 0)
  rownames(cluster_weight) = cluster_weight$Row.names
  cluster_weight = cluster_weight[,-1]
  
  percentage = as.data.frame(cluster_weight %>% group_by(y) %>% summarise(across(everything(), sum)))
  #colnames(percentage) = c( "Cluster",'Immune cell',"CAF","TAM","HPC-like cell","Malignant cell")
  colnames(percentage)[1] = c( "Cluster")
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
  #level_names <- c("Malignant cell","Immune cell","CAF","TAM","HPC-like cell")
  level_names <- sort(levels(percentage_long$variable))
  percentage_long$variable <- factor(percentage_long$variable, levels=level_names)
  
  percentage_long$Cluster <- factor(percentage_long$Cluster, levels = 1:12)
  
  datList[[slice_set]] <- percentage_long
  if (slice_set == 1){
    figure_list[[slice_set]] <- bar_plot(percentage_long,geom_bar_position=geom_bar_position, legend_position='right',
                                         color_pal =cols_cluster )
  }else{
    figure_list[[slice_set]] <- bar_plot(percentage_long,geom_bar_position=geom_bar_position, legend_position='none',
                                         color_pal = cols_cluster)
  }
}

if(need_scale){
  ggsave(file=paste0('result/Bulb_percentage_scaled_bulb',slice_set_sequence,'.png'), plot = figure_list[[slice_set]], 
         width = 10, height = 8, units = "in", bg = 'white', dpi = 500,limitsize = F)
}else{
  ggsave(file=paste0('result/Bulb_percentage_unscaled_bulb',slice_set_sequence,'.png'), plot = figure_list[[slice_set]], 
         width = 10, height = 8, units = "in", bg = 'white', dpi = 500,limitsize = F)
}





##### 4e #####
cols_alpha2 <- c(cols_alpha, '#57A6D9')
main_order_names2 <- c(main_order_names)
main_order_names2[1] <- iDR_SC_newname

## Plot Macfadden's adjusted R2
load('MacR2Vec_Bulb16_merge70.rds')
p1 <- barPlot_real(MacR2Vec[main_order_names2], cols = cols_alpha2, ylabel = "McFadden_R2")+theme_classic() + 
  theme( axis.text.x = element_blank() ) + theme(text = element_text(size=20))
ggsave(file="./barPlot_MacFaddenR2_bulb16_merge.png", plot = p1,
       width =7, height =4, units = "in", dpi = 800)



load("metric_clusterMat_deconCelltype_mergesubtype_bulb16_mergeV2.rds")

p1 <- volinPlot_real(ariMat[,main_order_names2], cols = cols_alpha2)
p2 <- volinPlot_real(nmiMat[,main_order_names2], ylabel = "NMI",  cols = cols_alpha2)
ggsave(file="./volinPlot_ariMat_bulb16_mergeV2.png", plot = p1,
       width =7, height =4, units = "in", dpi = 800)
ggsave(file="/volinPlot_nmiMat_bulb16_mergeV2.png", plot = p2,
       width =7, height =4, units = "in", dpi = 800)

##### 4f #####

#### Plot tractory on tSNE
load("pseudotime_allspots_bulb16_merge.rds")
load("tsne2_allMethods_Bulb16_rep2_merge70.rds")
tsne_idrsc <- tsne2List[[1]]

dat1 <- data.frame(PC1=tsne_idrsc[,1], PC2=tsne_idrsc[,2], pseudotime=pseudotime_allspots,
                   scaled_pseudotime = range01(pseudotime_allspots, na.rm=T))
library(ggplot2)
med <- quantile(dat1$scaled_pseudotime, 0.3, na.rm=T)
p1 <- ggplot(subset(dat1, !is.na(pseudotime)), aes(x=PC1, y=PC2, color=scaled_pseudotime)) + 
  geom_point( alpha=1) +
  scale_color_gradientn(colours = c( "#FFBB78", "#96B5E9", "#F54902"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size=14,color="black"),
        axis.text.y=element_text(size=14, color="black"),
        axis.title.x = element_text(size=18, color='black'),
        axis.title.y = element_text(size=18, color='black'),
        title= element_text(size=20, color='blue'),
        legend.text=element_text(size=14),
        legend.position = 'none',
        panel.background= element_rect(fill = 'white', colour = 'black')) +
  cowplot::theme_cowplot()+ xlab('tSNE 1') + ylab("tSNE 2") + theme(legend.position = 'none')
p1

ggsave(file=paste0("Bulb16_merge70_TSNEpsedotimeV2.png"), plot = p1,
       width = 6.5, height =5, units = "in", dpi = 400)


load("pseudotime_tmpList_bulb16_merge70.rds")
load('posList16_before_merge70_Bulb16rep2.rds')
posList <- posList16_before_merge70
sum(sapply(posList, nrow))



idx_outList <- list()
for(r in 3:16){
  idx_outList[[r]] <-  c(order(posList[[r]][,2], decreasing = T)[1:10],
                         order(posList[[r]][,1])[1:10])
  
}
for(r in 1:2){
  idx_outList[[r]] <- 25
}
pt_size = 0.1
pList_pseudo <- list()
datList <- list()
for(r in 1:16){ ## each sample
  # r <- 1
  
  
  message("r = ", r)
  pseudotime_tmp <- pseudotime_tmpList[[r]]  [-idx_outList[[r]]]
  indx <- 1:length(pseudotime_tmp)#which(!is.na(pseudotime_tmp))
  pseudotime_tmp <- pseudotime_tmp[indx]
  pos_tmp <- posList[[r]][-idx_outList[[r]],][indx,]
  
  dat1 <- data.frame(row=pos_tmp[,1], col=pos_tmp[,2], pseudotime=pseudotime_tmp,
                     scaled_pseudotime = range01(pseudotime_tmp, na.rm=T))
  dat1$Sample <- r
  med <- quantile(dat1$scaled_pseudotime, 0.6, na.rm=T)
  datList[[r]] <- dat1[,-3]
  ptmp <- ggplot(dat1, aes(x=row, y=col, color=scaled_pseudotime)) + geom_point(size=pt_size, alpha=1) +
    scale_color_gradientn(colours = c( "#FFBB78", "#96B5E9", "#F54902"))+
    mytheme_graybox()  + theme(legend.position = 'none') 
  pList_pseudo[[r]]<- ptmp
  
  
}
library(cowplot)
pList_pseudo[[1]]
p12 <- plot_grid(plotlist =pList_pseudo[1:8], nrow=2, ncol=4, byrow=T)
ggsave(file=paste0("Bulb16_1to8_merge70_psedotime_spatialheatmapV2.png"), plot = p12,
       width = 12, height =6, units = "in", dpi = 400)


# p12 <- plot_grid(plotlist =pList_pseudo[9:16], nrow=2, ncol=4, byrow=T)
# ggsave(file=paste0("Bulb16_9to16_merge70_psedotime_spatialheatmapV2.png"), plot = p12,
#        width = 12, height =6, units = "in", dpi = 400)
# get legend
ggplot(dat1, aes(x=row, y=col, color=scaled_pseudotime)) + geom_point(size=pt_size) +
  scale_color_gradientn(colours = c( "#FFBB78", "#96B5E9", "#F54902"))+ theme_bulb16()

