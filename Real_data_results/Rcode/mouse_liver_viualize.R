### dir_current <- "F:/Research paper/IntegrateDRcluster/AnalysisCode/PRECAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/utility_funcs.R"))
setwd(paste0(dir_current, "Real_data_results/dataFiles/Liver8/") )

##### 3a #####

load("metric_cluster_annotate_Liver8V2.rds")
apply(ariMat, 2, median)
cols_alpha2 <- c(cols_alpha, '#57A6D9')
main_order_names2 <- c(main_order_names)
main_order_names2[1] <- iDR_SC_newname

p1 <- volinPlot_real(ariMat[,main_order_names2], cols = cols_alpha2)

ggsave(file="./output_figs/volinPlot_ariMat_mouseLiver8.png", plot = p1,
       width =5, height =5, units = "in", dpi = 800)

p1 <- barPlot_real(ariVec[main_order_names2], cols = cols_alpha2)+theme_classic() + 
  theme( axis.text.x = element_blank() ) + theme(text = element_text(size=20))
ggsave(file="./output_figs/barPlot_ariVec_mouseLiver8.png", plot = p1,
       width =5, height =5, units = "in", dpi = 800)

load("MethodR4_batchRemovalMat_Liver8V2.rds")
load("MethodPy5_batchRemovalMat_Liver8.rds")
clisiMat <- cbind(clisiMat, clisiMat_py)
ilisiMat <- cbind(ilisiMat, ilisiMat_py)

p1 <- volinPlot_real(clisiMat[,main_order_names2], cols = cols_alpha2, ylabel='cLISI') + ylim(c(1,3))
ggsave(file="./output_figs/volinPlot_cLISI_mouseLiver8.png", plot = p1,
       width =5, height =5, units = "in", dpi = 800)
p2 <- volinPlot_real(ilisiMat[,main_order_names2], ylabel = "iLISI",  cols = cols_alpha2)
ggsave(file="./output_figs/volinPlot_iLISI_mouseLiver8.png", plot = p2,
       width =5, height =5, units = "in", dpi = 800)

##### 3b #####
load("hZtsneList_mouseLiver8.rds")

library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_pal <- c(pal1(max_n), pal2(max_n2)[c(1,3,8,12)]);
idx_used <- c(6,5, 1, 4, 2 , 8,7)
cols_cluster2 <- cols_pal[idx_used]

## select colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols_sample <- gg_color_hue(8)


pt_size_sample <- 0.5; pt_alpha_sample <- 0.8
pt_size_cluster <- 0.5; pt_alpha_cluster <- 1
base_axis_size <- 20

library(ggplot2)
str(hZtsneList)
pList_sample <- list(); pList_cluster <- list()
for(j in 1:length(hZtsneList)){
  message('j = ', j)
  
  
    
    pList_sample[[j]] <- plot_scatter(hZtsneList[[j]], meta_data, label_name = 'sample',border_col="gray10", 
                                      base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
                                      point_alpha = pt_alpha_sample,no_guides=T) + labs(x=NULL, y=NULL)
    pList_cluster[[j]] <- plot_scatter(hZtsneList[[j]], meta_data, label_name = 'cluster',
                                       base_size=20,border_col="gray10",
                                       palette_use=cols_cluster2,
                                       point_size = pt_size_cluster, 
                                       point_alpha =pt_alpha_cluster,no_guides=T) + labs(x=NULL, y=NULL)
    
  
  
  
  
}
names(pList_sample) <-names(pList_cluster) <-  names(hZtsneList)

library(cowplot)
#pList_sample[[1]]
subnames <- c("iDR-SC" ,"Seurat V3", "fastMNN","MEFISTO", "Uncorrected")
p12 <- plot_grid(plotlist =pList_sample[subnames], nrow=1, ncol=5)
ggsave(file="./output_figs/Liver8_tsne_sample_heatmap.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 200)
p1 <- pList_cluster[[1]]
pList_cluster[[1]] <- Seurat::LabelClusters(p1, id='cluster', repel = T, 
                                            size = 9, color = 'white', 
                                            box = T, position = "nearest")
p12 <- plot_grid(plotlist =pList_cluster[subnames], nrow=1, ncol=5)
ggsave(file="./output_figs/Liver8_tsne_cluster_heatmap.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 200)

##### 3c #####
library(Seurat)
load("seu_DEG_heatmap_combinedSample.rds")
seu <- seuAll
seu@assays$RNA@var.features <- row.names(seu)
seu <- ScaleData(seu)
cols_cluster <- c("#98df8a", "#2ca02c", "#1f77b4", "#ffbb78", "#aec7e8", "#ff9896", "#d62728")
color_id <- as.numeric(levels(Idents(seu)))

genes_use <- c("Cyp2e1", "Oat", "Cyp2c37", "Gulo", "Glul", "Slc1a2", 'Cyp2d9', "Gstm3", 
               "Malat1", "Cox1","Hamp2", "Cyp3a44", "Gsn", "Dpt","Vim", "Tagln",  
               "Cyp2f2", "Sds", "Hal", "Ctsc", "Aldh1b1", "Hsd17b13", "Spp1")


p1 <- doHeatmap(seu, features = genes_use, cell_label= "Domain",
                grp_label = F, grp_color = cols_cluster[color_id],
                #disp.max = 2.1,
                pt_size=6,slot = 'scale.data') + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=12, face= "italic", family='serif'))
ggsave(paste0("./output_figs/Liver8All_usedGenesDEGs_heatmapV2.png"), plot = p1, width = 8, height = 6, units = "in", dpi = 1000)

##### 3d #####
load('MacR2Vec_reference1_mouseLiver.rds')
names(MacR2Vec) <- main_order_names2
p1 <- barPlot_real(MacR2Vec[main_order_names2], cols = cols_alpha2, ylabel = "McFadden's R^2")+theme_classic() + 
  theme( axis.text.x = element_blank() ) + theme(text = element_text(size=20))
p1


##### 3e #####
load("metadata_filter2List_mouseLiverST.rds")
load("idrsc_cluster7_mouseLiver8.rds")
posList <- lapply(metadata_filter2List, function(x){
  y <- cbind(x$row, x$col)
  row.names(y) <- row.names(x)
  y
} )
color_pal = list("#CB181D","#EF3B2C","#FB6A4A","#FC9272","#FCBBA1","#1f77b4","#ff7f0e","#2ca02c","#c5b0d5")

# sample 1: slice 1 to 3
dat_List_all <- list()
figure_list = list()
for (slice_set in 1:3){
  # slice_set <- 1
  st_coord = posList[[slice_set]]
  
  norm_weights = read.csv(paste('./deconvolution_data/output_weights_mouseLiver8_st',slice_set,'.csv',sep=""),row.names=1)
  
  cluster_id = cbind(posList[[slice_set]],clusterList[[slice_set]])
  colnames(cluster_id) <- c("row", "col", "y")
  cluster_weight = merge(norm_weights,cluster_id,by = 0)
  rownames(cluster_weight) = cluster_weight$Row.names
  cluster_weight = cluster_weight[,-1]
  
  
  
  ct_set =  colnames(cluster_weight)[1:6]#[1:6] #"Hepatocyte_Alb.high" # "B.cell_Jchain.high" #
  head(cluster_weight)
  tmpList <- list()
  for(jj in 1: length(ct_set)){
    ## jj <- 1
    ct <- ct_set[jj]
    plot_val = as.data.frame(cluster_weight[,ct])
    rownames(plot_val) = rownames(cluster_weight)
    colnames(plot_val) = colnames(cluster_weight)[ct]
    barcodes = rownames(plot_val)
    my_table = as.data.frame(cluster_id[barcodes, c(1,2)] )
    my_table$Proportion = plot_val[barcodes,]
    my_table$cluster = cluster_weight[barcodes,'y']
    ylimit = c(0, 1)
    
    my_table$border_color = '#cccccc'
    my_table$stroke = 0
    
    select_cluster = c(1,2)
    for (i in select_cluster){
      my_table[my_table$cluster==i,'stroke'] = 0.45
    }
    
    for (i in 1:10){
      my_table[my_table$cluster==i,"border_color"] = color_pal[i]
    }
    
    table(my_table$cluster)
    my_table$Proportion <- range01(my_table$Proportion)
    if(ct ==  "Erythroid.cell_Hbb.bs.high"){
      med <- quantile(my_table$Proportion, 0.987)
    }else{
      med <- quantile(my_table$Proportion, 0.94)
    }
    
    my_table$Proportion[ my_table$Proportion > med] <- med
    my_table$Proportion <- range01(my_table$Proportion)
    
    tmpList[[jj]] <- my_table
    plot <- ggplot(my_table, aes(x = row, y = col)) + 
      geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 5,pch=21) +
      scale_fill_gradientn(colors = c("#361A95","white", "#D62728")) +  # ,limits = ylimit,#FBA536
      # scale_colour_gradient2(
      #   low = "#361A95",
      #   mid = "white",
      #   high = "#FBA536", midpoint = med)+
      scale_shape_identity() + scale_size_identity() + theme_classic() +
      mytheme_graybox(base_size = 28) + scale_y_reverse() + theme(legend.position = 'none')
    
    figure_list[[((slice_set-1) * length(ct_set) + jj)]] <- plot #+ coord_fixed() # + xlim(xlim) + ylim(ylim)
    
  }
  dat_List_all[[slice_set]] <- tmpList
}



p12 <- cowplot::plot_grid(plotlist=figure_list, nrow=3, ncol=6, byrow=T)
ggsave(file='output_figs/select_type_Sample1.png', plot = p12, 
       width = 14, height = 6, units = "in", bg = 'white', dpi = 50,limitsize = F)

##### 3f #####

load("pseudotime_allspots_mouseLiver8V2.rds")
load("hZtsneList_mouseLiver8.rds")
pseudo.slingshot <- pseudo.slingshot2
str(pseudo.slingshot2)
tsne_idrsc <- hZtsneList[[1]]

str(tsne_idrsc)
dat1 <- data.frame(PC1=tsne_idrsc[,1], PC2=tsne_idrsc[,2], pseudotime=pseudo.slingshot,
                   scaled_pseudotime = range01(pseudo.slingshot, na.rm=T))
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
  cowplot::theme_cowplot()+ xlab('tSNE 1') + ylab("tSNE 2")# + theme(legend.position = 'none')
p1

##### 3g #####
library(SingleCellExperiment)
load("Liver8_sce_changedGenes_along_pseudotime.rds")
load('Liver8_pseudo_TSCANtest.rds')
load("pseudotime_allspots_mouseLiver8V2.rds")
pt <- pseudo.slingshot2
pseudo[order(pseudo$p.value),]
#sum(abs(pseudo$logFC) > 0.01)
sum(pseudo$FDR<0.05)
sorted <- pseudo[order(pseudo$p.value),]
up.left <- sorted # [sorted$logFC < 0,]
head(up.left, 10)
best <- head(row.names(up.left), 20)

rowData(sce.liver)$SYMBOL <- row.names(sce.liver)
sce.liver$Pseudotime <- pt

domain_colors <- c("#98df8a", "#2ca02c", "#1f77b4", "#ffbb78", "#aec7e8", "#ff9896", "#d62728")
domain_names_allspots <- cl_idrsc_cellRegion

## orered top 30 genes
features_reorder <- c("Glul", "Oat","Cyp2e1","Slc1a2","Gulo","Lect2" ,  "Cyp2e1",
                      "Cyp2c37","Cyp2c69", "Cyp3a11","Rnase4", "Aldh1a1","Cyp2c29",
                      ## lower
                      "Cyp2f2", "Hal" ,"Sds", "Hsd17b13", "Ctsc", "A1bg","Rida")

p1 <- scater::plotHeatmap(sce.liver, order_columns_by="Pseudotime", 
                          colour_columns_by= "domain", features= features_heat,
                          center=TRUE, swap_rownames="SYMBOL", cluster_rows=F)


