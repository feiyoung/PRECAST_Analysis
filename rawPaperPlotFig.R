
dir_current <- "E:\\Research paper\\IntegrateDRcluster\\PRECAST_Analysis\\"

source(paste0(dir_current, "//utility_funcs.R"))
##### 2a #####
setwd(paste0(dir_current, "dataFiles\\brain12\\") )

library(cowplot)
load("Brain12_RGB_hZ_umapList_allsample.rds")
load("brain12_posList.rds")
load("RGB_umap_tsne_paste_brain12.rds")
str(hZ_umap3List_allsample)
## Add the results of PASTE
hZ_umap3List_allsample[["PASTE"]] <- hZ_umap3_paste

### RGB plot for each sample based Joint UMAP
indexList <- get_indexList(posList)
hZ_umap3List_allsample[2:3] <- hZ_umap3List_allsample[3:2]
names(hZ_umap3List_allsample)[2:3] <- names(hZ_umap3List_allsample)[3:2]
### Plot 151674, r=10
pList_umap_RGB_r10 <- list()
for(i in 1:9){ ## each method
  message("i = ", i)
  r <- 10
  ptmp <- plot_RGB(posList[[r]], hZ_umap3List_allsample[[i]][indexList[[r]], ], pointsize = 0.4) + mytheme_graybox()
  ptmp <- rotate_90_clockwise(ptmp) 
  pList_umap_RGB_r10[[i]] <- ptmp
}
names(pList_umap_RGB_r10) <- names(hZ_umap3List_allsample)
subnames <- c("iDR-SC",  "Seurat V3",  "Harmony",    "fastMNN")
p12 <- plot_grid(plotlist =pList_umap_RGB_r10[subnames], nrow=1, ncol=4)



load("Brain12_clusterPyMethod4_res.rds")
clustM1 <- clusterMat
load("Brain12_clusterRMethod_res.rds")

str(clusterMat)
head(clusterMat)
clusterMat <- cbind(clusterMat[,1:5], clustM1)
head(clusterMat)
colnames(clusterMat)[2] <- "Seurat V3"
load("Brain12_metric_cluster10_idrsc.rds")
clusterMat <- cbind("iDR-SC"=unlist(cluster_idrsc), clusterMat)

load("cluster_allsample_paste_brain12.rds")
str(yMat_paste)
clusterMat <- cbind(clusterMat, "PASTE-SC-MEB"=yMat_paste[,1])

MethodNames <- c("iDR-SC", "Seurat V3" ,"Harm-SC-MEB",
                 "fastMNN-SC-MEB",  "Scanorama-SC-MEB",
                 "scGen-SC-MEB", "scVI-SC-MEB", 'MEFISTO-SC-MEB', "PASTE-SC-MEB")
clusterMat_sub <- clusterMat[, MethodNames]
dim(clusterMat_sub)
load(file='brain12_yList.rds')
yList <- lapply(yList, as.numeric)
for(r in 1:12){
  y <- yList[[r]]
  y[is.na(y)] <- 7
  yList[[r]] <- factor(y)
}
clusterMat_sub <- cbind(clusterMat_sub, Groundtruth = unlist(yList))


load("brain12_posList.rds")

indexList <- get_indexList(posList)
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names(palettes)
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_cluster <- c(tableau_color_pal()(10), pal1(7)[-1])
#simutool::colorbar_adj_transparent(cols_cluster[1:10], alpha = 1)

base_cols_cluster <- 1:7
base_cluster <- clusterMat_sub[,10]
ClusterMat <- clusterMat_sub[,1:9]


colorList <- c(align_colors(ClusterMat, base_cluster, cols_cluster), 
               list(cols_cluster[1:7]))
cols_cluster_idrsc <- colorList[[1]]
cols_cluster_idrsc[c(9, 10)] <- cols_cluster_idrsc[c(10, 9)]
colorList[[1]] <- cols_cluster_idrsc  
#simutool::colorbar_adj_transparent(cols_cluster_idrsc, alpha = 1)



### plot brain10
pList_assign_r10 <- list()
apply(clusterMat_sub, 2, max)
sapply(colorList, length)
for(j in 1:10){ ## each method
  # j <- 1
  pList_eachmethod <- list()
  for(r in 10){ ## each sample
    message("r = ", r)
    cluster_tmp <- clusterMat_sub[indexList[[r]],j]
    p_tmp <- plot_scatter(posList[[r]], meta_data = data.frame(cluster=factor(cluster_tmp)),
                          label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                          point_size = 1,palette_use = colorList[[j]], point_alpha = 0.6) + 
      mytheme_graybox()
    
    
    pList_assign_r10[[j]]<- rotate_90_clockwise(p_tmp)
  }
}
names(pList_assign_r10) <- colnames(clusterMat_sub)

subnames <- c("iDR-SC" ,"Seurat V3",   "Harm-SC-MEB" ,"fastMNN-SC-MEB")
library(cowplot)
p12 <- plot_grid(plotlist =pList_assign_r10[subnames], nrow=1, ncol=4)



p_annotation <- pList_assign_r10[["Groundtruth"]]
ggsave(file="brain12_assign_heatmap_annotation_r10.png", plot = p_annotation,
       width = 3, height =2, units = "in", dpi = 500)


##### 2b #####

load("brain12_tsne2_allMethods.rds")
load("Brain12_metric_cluster10_idrsc.rds")
meta_data$cluster_idrsc <- factor(unlist(cluster_idrsc), levels=1:10)
load("tsne2_paste_brain12.rds")
tsne2List[["PASTE"]] <- hZ_tsne2_paste

library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names(palettes)
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_pal <- c(pal1(max_n), pal2(max_n2)[c(1,3,8,12)]);

# plot(1:3, col=scales::hue_pal()(3), lwd=5)
library(RColorBrewer)
display.brewer.all() 
cols_sample <- c(
  brewer.pal(8,"Reds")[seq(2,8, by=2)],
  brewer.pal(8,"Greens")[seq(2,8, by=2)],
  brewer.pal(8,"Blues")[seq(2,8, by=2)]
)
base_cols <- gg_color_hue(3)
cols_sample[4] <- base_cols[1]
cols_sample[8] <- base_cols[2]
cols_sample[12] <- base_cols[3]

# global settings: 
pt_size_sample <- 0.2; pt_alpha_sample <- 0.7
pt_size_cluster <- 0.6; pt_alpha_cluster <- 0.5
base_axis_size <- 20
cols_cluster <- c( "#E15759CC", "#59A14FCC", "#FF9DA7CC", "#9C755FCC",
                   "#B07AA1CC", "#4E79A7CC", "#EDC948CC", "#BAB0ACCC",
                   "#ff7f0e", "#aec7e8CC")
library(ggplot2)

pList_sample <- list(); pList_cluster <- list()
for(j in 1:length(tsne2List)){
  message('j = ', j)
  pList_sample[[j]] <- plot_scatter(tsne2List[[j]], meta_data, label_name = 'sample',border_col="gray10",  
                                    base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
                                    point_alpha = pt_alpha_sample,no_guides=T) + labs(x=NULL, y=NULL)
  pList_cluster[[j]] <- plot_scatter(tsne2List[[j]], meta_data, label_name = 'cluster_idrsc', base_size=20,
                                     palette_use=cols_cluster, point_size = pt_size_cluster, border_col="gray10", 
                                     point_alpha =pt_alpha_cluster,no_guides=T) + labs(x=NULL, y=NULL)
  
}
names(pList_sample) <-names(pList_cluster) <-  names(tsne2List)
#pList_cluster[[1]]
library(cowplot)
subnames <- c("iDR-SC" ,"Seurat V3", "Harmony","fastMNN", "Uncorrected")

p12 <- plot_grid(plotlist =pList_sample[subnames], nrow=1, ncol=5)

ggsave(file="brain12_tsne_sample_heatmap.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 500)

p1 <- pList_cluster[[1]]
pList_cluster[[1]] <- Seurat::  LabelClusters(p1, id='cluster_idrsc', repel = T, 
                                              size = 11, color = 'white', 
                                              box = T, position = "nearest")

p12 <- plot_grid(plotlist =pList_cluster[subnames], nrow=1, ncol=5)
ggsave(file="brain12_tsne_cluster_heatmapV2.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 500)

## Get legends
plot_scatter(tsne2List[[1]], meta_data, label_name = 'sample',border_col="gray10",  
             base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
             point_alpha = pt_alpha_sample,no_guides=F)
plot_scatter(tsne2List[[1]], meta_data, label_name = 'cluster_idrsc', base_size=20,
             palette_use=cols_cluster_idrsc, point_size = pt_size_cluster, border_col="gray10", 
             point_alpha =pt_alpha_cluster,no_guides=F)

##### 2c #####
## boxplot for each sample
load('Brain12_eachsample_cluster_allMethods_res.rds')
MethodNames <- c("iDR-SC", "Seurat V3" ,"Harm-SC-MEB",
                 "fastMNN-SC-MEB",  "Scanorama-SC-MEB",
                 "scGen-SC-MEB", "scVI-SC-MEB", 'MEFISTO-SC-MEB')
ariMat_sub <- ariMat[,MethodNames]
apply(ariMat_sub, 2, median)
colnames(ariMat_sub) <- main_order_names[-9]

## Add PASTE

load("./paste_ARI_NMI_allsample_brain12.rds")
ari_paste <- cluster_meaMat_paste_allsample[,1]
colnames(ariMat_sub)[1] <- iDR_SC_newname
p1 <- volinPlot_real(cbind(ariMat_sub, "PASTE"=ari_paste), cols=cols_alpha2 ) + theme(legend.position = 'bottom')


##### 2d #####
load("hZ_idrsc_brain12.rds")
load("Brain12_metric_cluster10_idrsc.rds")
cluster_idrsc_vec <- unlist(cluster_idrsc)

sublist <- subsample(hZ_idrsc, cluster_idrsc_vec, 0.010)
order_idx <- order(sublist[[2]])
cluster_orderd <- sublist[[2]][order_idx]
corMat <- cor(t(sublist[[1]][order_idx, ]))
row.names(corMat) <- paste0('sp', 1:nrow(corMat))
colnames(corMat) <- paste0('spot', 1:nrow(corMat))
corMat[1:4,1:5]
p1 <- doHeatmap.matrix(corMat, cluster_orderd, legend_title = 'Domain', grp_color=cols_cluster)

ggsave(file="hZ_idrsc_cor_heatmap_brain12.png", plot = p1,
       width =10, height =8, units = "in", dpi = 100)

##### 2e #####
gene_eachdomain <- c("HOPX", "PCP4", "Enc1", "CNP", "MBP", "NEFL")
names(gene_eachdomain) <- paste0("Domain", c(1,2,4, 5, 7, 8))
gene_eachdomain <- firstup(gene_eachdomain)
cols <- c('#3AB370',"#FD1593") # 
# simutool::colorbar_adj_transparent(cols,alpha = 1)

## plot for sample 10
load("brain_seu_sample10.rds")
title_size <- 14
pt_size <- 1
quantVec <- c(0.92, 0.95, 0.92, 0.93, 0.9, 0.96)
p2List <- list()
for(j in 1: length(gene_eachdomain)){
  # j <- 5
  message("j = ", j)
  p1 <- featurePlot(seu, feature = gene_eachdomain[j], pt_size=pt_size, quant=quantVec[j],
                    cols = cols ) + ggtitle("")
  
  
  p2List[[j]] <- rotate_90_clockwise(p1)
}
p12 <- cowplot::plot_grid(plotlist = p2List, nrow=1, ncol=6)
#p12
ggsave(file="spaHeatmap_DEGs_brain_sample10V3.png", plot = p12,
       width = 13, height =2.5, units = "in", dpi = 500)


p1 <- featurePlot(seu, feature = gene_eachdomain[1], pt_size=pt_size, quant = 0.97,
            cols = cols ) + theme(legend.position = 'right')
rotate_90_clockwise(p1)


##### 2f #####
corrSig <- read.csv(file='corrTestMat_sig.csv', row.names = 1)
head(corrSig)
cutoff <- 0.001

row.names(corrSig) <- firstup(row.names(corrSig))
corrSig_new <- corrSig[corrSig[,4]<cutoff,]

gene_assoc_pseudotime <- row.names(corrSig_new)[order(abs(corrSig_new[,1]), decreasing = T)][1:6]

quantVec <- c(0.95, 0.95, 0.95, 0.95, 0.9, 0.9)
p2List <- list()
for(j in 1: length(gene_assoc_pseudotime)){
  # j <- 1

  p1 <- featurePlot_pseudoassoc(seu, feature = gene_assoc_pseudotime[j], pt_size=1, quant=quantVec[j]
  ) + 
    ggtitle("")
  
  p2List[[j]] <- rotate_90_clockwise(p1)
}
p12 <- cowplot::plot_grid(plotlist = p2List, nrow=1, ncol=6)
p12
ggsave(file="spaHeatmap_pseudotimeGenes_brain_sample10.png", plot = p12,
       width = 13, height =2.5, units = "in", dpi = 500)



##### 2g #####
library(org.Hs.eg.db)
library(limma)
Gid_tmp <- mapIds(org.Hs.eg.db, toupper(row.names(corrSig)), column='ENTREZID',
                  keytype ='SYMBOL')
goa.DP.down <- goana(Gid_tmp, species = "Hs")
head(goa.DP.down)
pv005 <- goa.DP.down$P.DE[goa.DP.down$P.DE<1]
padj <- p.adjust(pv005, method="fdr")
sum(padj<0.05)
n_BP <- sum(goa.DP.down$Ont == 'BP')
n_CC <- sum(goa.DP.down$Ont == 'CC')
n_MF <- sum(goa.DP.down$Ont == 'MF')
set.seed(1)
perb_BP <- sample(n_BP)
perb_CC <- sample(n_CC) + n_BP
perb_MF <- sample(n_MF) + n_BP + n_CC
df1 <- rbind(goa.DP.down[perb_BP,], goa.DP.down[perb_CC,], goa.DP.down[perb_MF,])
df1$nlog10P <- -log10(df1$P.DE)
head(df1, 5)
id_label_bp <- order(subset(df1, Ont=='BP')[,'nlog10P'], decreasing = T)[1:3]
id_label_cc <- order(subset(df1, Ont=='CC')[,'nlog10P'], decreasing = T)[1:2]
id_label_mf <- order(subset(df1, Ont=='MF')[,'nlog10P'], decreasing = T)[1:2]

df1$ID <- 1:nrow(df1)
library(ggplot2)
p = ggplot(df1,aes(ID, nlog10P, color=Ont))+
  geom_point(aes(size=DE, color=Ont), alpha=0.5) +labs(color='Category', 
                                                       x="Pathway ID",y=expression(-log[10](P)),title=NULL)

cols <- c("steelblue3", "goldenrod", "brown3")
label <- rbind(subset(df1, Ont=='BP')[id_label_bp[c(1,2,3)],],
               subset(df1, Ont=='CC' )[id_label_cc,],
               subset(df1, Ont=='MF' )[id_label_mf,])
p1 <- p + cowplot::theme_cowplot() +
  ggrepel::geom_text_repel(data = label, aes(label = Term), size=6) +
  theme(axis.text.x=element_blank(),axis.title.x = element_blank() ) +
  scale_color_manual(values = cols) + 
  geom_hline(yintercept=-log10(0.05), linetype="dotted", color='black') + 
  guides(color=guide_legend("Category", override.aes = list(size = 4)))+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16),
        axis.title.y = element_text(size=18),
        legend.box = 'vertical', legend.box.just = 'top', legend.position = 'right')
p1
ggsave(file="GO_Enrichment_pseudogenes_brain12.png", plot =p1,
       width = 13, height =2.5, units = "in", dpi = 500)


##### 2f #####
load('profiler_termList_brain12.rds')
library(gprofiler2)
library(ggplot2)
source_set <- c('KEGG')
pathway_remove <- "KEGG:05171"
pList <- list()

r <- 10
cat("r = ", r, '\n')
gostres2 <- termList[[r]]
gostres2$result <- subset(gostres2$result, source %in% source_set)
gostres2$meta$query_metadata$sources <- source_set
p <- gostplot(gostres2, capped = FALSE, interactive = F)  + theme(text=element_text(size=20))
label_terms <- get_top_pathway1(gostres2$result, ntop=6, source_set = source_set)$term_id
label_terms <- label_terms[label_terms!= pathway_remove]
p1 <- publish_gostplot(p, highlight_terms = label_terms,
                       filename=paste0("sample10_KEGG_SVGs_pathways.pdf"), 
                       width = 9, height =6)



##### 3b #####
rm(list=ls())

setwd(paste0(dir_current, "dataFiles\\HCC4\\") )
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
for(r in 1:4){## each sample
  # i <- 1
  
  pList_eachmethod <- list()
  for(i in 1){  ## each method
    # r <- 1
    message("r = ", r)
    pt_size = 0.5
    if(r == 1) pt_size= 0.3
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



##### 3c #####

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


##### 3d & 3e #####
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

setwd("D:\\LearnFiles\\Research paper\\IntegrateDRcluster\\RealData\\Rcode\\HCC4Data")
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
    geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 5,pch=21,colour = my_table$border_color) +
    scale_fill_gradientn(colors = my_pal,limits = ylimit) + 
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
    geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 3,pch=21,colour = my_table$border_color) +
    scale_fill_gradientn(colors = my_pal,limits = ylimit) + 
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
    geom_point(aes(fill = Proportion),size = 2.3, stroke = my_table$stroke,alpha = 3.5,pch=21,colour = my_table$border_color) +
    scale_fill_gradientn(colors = my_pal,limits = ylimit) + 
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



##### 4b #####

setwd(paste0(dir_current, "dataFiles\\Bulb16\\") )
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
    ptmp <- plot_scatter(pos_tmp[idx, ], 
                         meta_data = data.frame(cluster=clusterlabel),
                         label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                         point_size = pt_size,palette_use = cols_cluster[unique_sort], point_alpha = point_alpha) + 
      theme_bulb16()
  }else{
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
for(r in 1:16){ ## each sample
  # r <- 1
  
  
  message("r = ", r)
  pseudotime_tmp <- pseudotime_tmpList[[r]]  [-idx_outList[[r]]]
  indx <- 1:length(pseudotime_tmp)#which(!is.na(pseudotime_tmp))
  pseudotime_tmp <- pseudotime_tmp[indx]
  pos_tmp <- posList[[r]][-idx_outList[[r]],][indx,]
  
  dat1 <- data.frame(row=pos_tmp[,1], col=pos_tmp[,2], pseudotime=pseudotime_tmp,
                     scaled_pseudotime = range01(pseudotime_tmp, na.rm=T))
  med <- quantile(dat1$scaled_pseudotime, 0.6, na.rm=T)
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

##### 5a #####

### 5a upper panel:

setwd(paste0(dir_current, "dataFiles\\Hippo2\\") )
load("hZ_umap3_tmpList_Method8_Hip2_merge70.rds")
str(hZ_umap3_tmpList)
## Add PASTE
load("paste_umap3_tsne3_tmpList_Hip2_merge70.rds")
for(r in 1:2){
  hZ_umap3_tmpList[[r]][["PASTE"]] <- umap3_paste_tmpList[[r]]
}

load('posList_raw_Hip2.rds')
idx_outList <- list()
idx_outList[[1]] <-  5
idx_outList[[2]] <- order(posList_raw[[2]][,1])[1:4]

nMethod <- 9
pList_umap_RGB <- list()
for(r in 1:2){ ## each sample
  # r <- 1
  if(r ==1) pt_size = 0.25
  if(r == 2) pt_size = 0.1
  message("r = ", r)
  for(i in 1:nMethod){ ## each method
    # i <- 1
    if(i <= nMethod){
      uMat <- hZ_umap3_tmpList[[r]][[i]][-idx_outList[[r]],]
      indx <- complete.cases(uMat)
      uMat <- uMat[indx,]
      ptmp <- plot_RGB(coordinate_rotate(posList_raw[[r]][-idx_outList[[r]],][indx,], rotate_angles[r]), 
                       uMat, 
                       pointsize = pt_size) + theme_hip2()+   scale_x_reverse()
      
      
      pList_umap_RGB[[(r-1)*nMethod+ i]]<- ptmp
    }
    
  }
  
}
library(cowplot)
# p12 <- plot_grid(plotlist =pList_umap_RGB, nrow=2, ncol=nMethod, byrow=T)
# ggsave(file="Hip2_merge70_umapRGB.png", plot = p12,
#        width = 14, height =4, units = "in", dpi = 200)



## for two samples of PRECAST
jmethod <- 1
method_index <- c(jmethod, 9+jmethod)
p12 <- plot_grid(plotlist = pList_umap_RGB[method_index], nrow=1, ncol=2)
ggsave(file=paste0("Hip2_merge70_umap_RGB_iDRSC.png"), plot = p12,
       width = 6, height =3, units = "in", dpi = 400)

### cluster assignment spatial heatmap
load("clusterAll_tmpList_Method8_Hip2_merge70.rds")
str(clusterAll_tmpList[[1]])
load("clusterAll_paste_tmpList_Hip2_merge70.rds")
for(r in 1:2){
  clusterAll_tmpList[[r]][["PASTE"]] <- clusterAll_paste_tmpList[[r]]
}

load('posList_raw_Hip2.rds')
idx_outList <- list()
idx_outList[[1]] <-  5
idx_outList[[2]] <- order(posList_raw[[2]][,1])[1:4]


##### domain color selection
library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names(palettes)
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")

cols_cluster <- c( "red2",'green3',"blue", c("#F4737A", "#E8818C", "#DC8F9E", "#D09DB0", "#C5ABC3"),
                   "#8FBC8F99", c('gold1',paste0("goldenrod",1:2), "goldenrod"),
                   c("darkolivegreen1", "darkolivegreen3"),c("darkorchid1", "darkorchid3"), pal1(max_n))
cols_cluster_low_reso <- cols_cluster

### color selection for each method mathched with iDR-SC
load("mergedSquare70_clusterList_chooseK17_Hip2.rds")
load("./clusterR3_SCMEB_methods7_Hip2_merge70.rds")

## change the order of clusters
rawID <- c(10, 11, 1, 7, 17, 8, 15, 13, 3, 9,  4,  16,2, 6, 5, 14, 12)
names(rawID) <- 1: 17
### Plot the domain heatmap of iDR-SC
pList <- list()
library(Seurat)
for(r in 1:2){
  #r <- 1
  i <- 1
  cluster_tmp <- clusterAll_tmpList[[r]][[i]][-idx_outList[[r]]]
  idx <- !is.na(cluster_tmp)
  cluster_tmp[idx] <- as.numeric(replace_ID(cluster_tmp[idx], rawID))
  unique_sort <- sort(unique(cluster_tmp))
  clusterlabel <- factor(cluster_tmp, levels = unique_sort)
  pos_tmp <- posList_raw[[r]][-idx_outList[[r]],]
  if(r==1) pt_size <- 0.8
  if(r==2) pt_size <- 0.5
  # -40.5 rotation and x reverse for r=2
  p1 <- plot_scatter(coordinate_rotate(pos_tmp[idx, ], rotate_angles[r]), 
                     meta_data = data.frame(cluster=clusterlabel[idx]),
                     label_name = 'cluster', xy_names = c("", ""), no_guides = T,
                     point_size = pt_size,palette_use = cols_cluster, point_alpha = 0.5) + 
    theme_hip2() + scale_x_reverse()
  pList[[r]] <- p1
  
}
library(cowplot)
p12 <- plot_grid(plotlist =pList, nrow=1, ncol=2)
ggsave(file="clusterAssign_iDRSC_Hip2_merge70V2.png", plot = p12,
       width = 6, height =3, units = "in", dpi = 400)


### 5a lower upper:


load('Hip2_raw_iDRSC_m24_cluster25.rds')
load("cluster_SCMEB_methodsR3_Hip2_RawRevo.rds")
load("posList_Hip2_RawRevo.rds")
cM1 <- cbind("iDR-SC"=unlist(cluster_idrsc_raw), clusterMat)
load("cluster_SCMEB_methods7_Hip2_RawRevo.rds")
clusterMat <- cbind(cM1[,1:4], clusterMat[,4:7])
apply(clusterMat, 2, max)
MethodNames <- c("iDR-SC", "Seurat V3" ,"Harm-SC-MEB",
                 "fastMNN-SC-MEB",  "Scanorama-SC-MEB",
                 "scGen-SC-MEB", "scVI-SC-MEB", 'MEFISTO-SC-MEB')
clusterMat_sub <- clusterMat[, MethodNames]
dim(clusterMat_sub)
indexList <- get_indexList(posList)
clusterMat <- clusterMat_sub
rawID <- c(25, 15, 17, 19, 4, 10, 13,  5, 6, 1, 21, 2, 14, 23,  22, 24, 3, 7, 8, 9, 11, 12, 16, 18, 20)
names(rawID) <- 1:25
clusterMat[,1] <- as.numeric(replace_ID(clusterMat[,1], rawID))

library(ggthemes)
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]
names(palettes)
pal1 <- tableau_color_pal("Classic 20")
max_n <- attr(pal1, "max_n")
pal2 <- tableau_color_pal("Classic Blue-Red 12")
max_n2 <- attr(pal2, "max_n")
cols_cluster = c( "#FD7446" ,"#709AE1", "#9EDAE5",
                  "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,
                  "#FF9896","#91D1C2", "#C7E9C0" ,
                  "#6B6ECF", "#7B4173",  pal1(max_n)) # "green4", 'blue',
cols_cluster[c(5, 15, 17, 25)] <- c( "darkseagreen","green3", "blue","red2")
cols_cluster[1:25] <- cols_cluster[1:25][rawID] 
cols_cluster[c(4:8)] <- c("#F4737A", "#E8818C", "#DC8F9E", "#D09DB0", "#C5ABC3")
cols_cluster[c(9:10)] <- c("#8FBC8F99", "#8FBC8F99") #c('tan1', RColorBrewer::brewer.pal(4, "YlOrRd")[2])
cols_cluster[c(11:19)] <- #RColorBrewer::brewer.pal(9, "Purples")
  c(c('gold1',paste0("goldenrod",1), "goldenrod"),#c("darkorchid1", "darkorchid3", "darkorchid4"),
    color2bar_gradient('steelblue1', 'steelblue3', 3), 
    color2bar_gradient('pink', 'pink3', 3))
cols_cluster[20] <- "darkolivegreen1"
cols_cluster[c(21:23)] <- c(color2bar_gradient("#8FBC8F99", 'green4', 3))#, #c("darkorchid1", "darkorchid3", "darkorchid4") #RColorBrewer::brewer.pal(6, "Greys")[c(3,5, 6)]
cols_cluster[c(24:25)] <- c("bisque2", 'bisque3')


base_cols_cluster <- 1:25
names(base_cols_cluster) <- cols_cluster[1:25]
base_cluster <- clusterMat[,1]
ClusterMat <- clusterMat[,2:8]



colorList <- c(list(cols_cluster[1:25]), 
               align_colors(ClusterMat, base_cluster, cols_cluster))

# simutool::colorbar_adj_transparent(cols_cluster[1:25], alpha = 1)
## find the outlinear point latter
idx_outList <- list()
idx_outList[[1]] <- 5
idx_outList[[2]] <- order(posList[[2]][,1])[1:2]

## ----------------------Plot PRECAST

remove_cluster <- 17:25
cluster_id_plot <- setdiff(1:25, remove_cluster)

cluster_idrsc <- clusterMat[,1]
cluster_idrsc_tmp <- cluster_idrsc
pList <- list()
for(r in 1:2){
  # r <- 2
  pt_size_set <- c(0.8, 0.5)
  cluster_tmp <- cluster_idrsc_tmp[indexList[[r]]][-idx_outList[[r]]]
  #cluster_tmp <- replace_ID(cluster_tmp, rawID)
  idx <- which(cluster_tmp %in% cluster_id_plot)
  # cluster_tmp[-idx] <- NA
  unique_sort <- sort(unique(cluster_tmp))
  clusterlabel <- factor(cluster_tmp)
  
  pos_tmp <- posList[[r]][-idx_outList[[r]],]
  
  dat1 <- data.frame(coordinate_rotate(pos_tmp, rotate_angles[r]))
  colnames(dat1) <- c("row", "col")
  ptmp <- ggplot(data=dat1, aes(x=row, y=col)) + geom_point(fill='gray', color='gray', size=0.05, alpha=0.5)+ 
    mytheme_graybox()+  scale_x_reverse()
  tmpdat <- cbind(coordinate_rotate(pos_tmp[idx, ], rotate_angles[r]),
                  as.numeric(clusterlabel[idx]))
  order_idx <- order(tmpdat[,3], decreasing = T)
  colnames(tmpdat) <- c("row1", "col1", "cluster")
  tmpdat <- as.data.frame(tmpdat)
  tmpdat$cluster <- factor(tmpdat$cluster)
  tmpp <- ptmp + geom_point(data= tmpdat[order_idx,], aes(x=row1, y=col1,color=cluster), size=pt_size_set[r]/4)+
    scale_color_manual(values=cols_cluster[cluster_id_plot])+theme_hip2() + theme(legend.position = 'none')
  
  pList[[r]] <- tmpp
}

p12 <- cowplot::plot_grid(plotlist=pList, nrow=1, ncol=2)
ggsave(file="Assign_iDRSC_Hip2.png", plot = p12,
       width = 6, height =3, units = "in", dpi = 400)




load("Hip2RawRevo_RGB_hZ_umapList_allsample.rds")
str(hZ_umap3List_allsample)
load("posList_Hip2_RawRevo.rds")
load("umap_tsne3_idrsc_mPotts24_Hip2.rds")
hZ_umap3List_allsample[[1]] <- umap3_idrsc

str(posList)
# change the location of Harmony and Seurat V3
indexList <- get_indexList(posList)
hZ_umap3List_allsample[2:3] <- hZ_umap3List_allsample[3:2]
names(hZ_umap3List_allsample)[2:3] <- names(hZ_umap3List_allsample)[3:2]

idx_outList <- list()
idx_outList[[1]] <-  5
idx_outList[[2]] <- order(posList[[2]][,1])[1:2]

## iDR-SC UMAP
pList_umap_RGB <- list()
for(r in 1:2){## each sample
  # r <- 1
  
  message("r = ", r)
  if(r == 1) pt_size <- 0.2
  if(r == 2) pt_size <- 0.1
  
  pt_size_set <- c(0.8, 0.5)
  cluster_tmp <- cluster_idrsc_tmp[indexList[[r]]][-idx_outList[[r]]]
  
  idx <- which(cluster_tmp %in% cluster_id_plot)
  # cluster_tmp[-idx] <- NA
  unique_sort <- sort(unique(cluster_tmp))
  clusterlabel <- factor(cluster_tmp)
  
  pos_tmp <- coordinate_rotate(posList[[r]][-idx_outList[[r]],], rotate_angles[r])
  
  tmpdat <- cbind(coordinate_rotate(pos_tmp[idx, ], rotate_angles[r]),
                  as.numeric(clusterlabel[idx]))
  order_idx <- order(tmpdat[,3], decreasing = T)
  
  
  idx_tmp <- indexList[[r]][-idx_outList[[r]]]
  ptmp <- plot_RGB(pos_tmp[idx,][order_idx,], hZ_umap3List_allsample[[1]][idx_tmp[idx], ][order_idx,],
                   pointsize = pt_size) +theme_hip2()+scale_x_reverse() + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  
  pList_umap_RGB[[r]]<- ptmp 
  
  
}

p12 <- plot_grid(plotlist = pList_umap_RGB, nrow=1, ncol=2)
ggsave(file=paste0("Hip2_RawRevo_umap_RGB_iDRSC.png"), plot = p12,
       width = 6, height =3, units = "in", dpi = 400)


##### 5b #####

### tSNE plots vs sample and clusters
load("hip2_merge70_tsne2_allMethods.rds")
str(tsne2List)
head(meta_data)
rawID <- c(10, 11, 1, 7, 17, 8, 15, 13, 3, 9,  4,  16,2, 6, 5, 14, 12)
names(rawID) <- 1: 17
meta_data$cell_type <- factor(replace_ID(meta_data$cell_type, rawID), levels=1:17)
load("tsne2_paste_hip2_merge70.rds")
tsne2List[["PASTE"]] <- hZ_tsne2_paste

cols_sample <- gg_color_hue(2)

cols_cluster2 <- cols_cluster_low_reso
cols_cluster2[2] <- 'green3'
colorbar_adj_transparent(cols_cluster2)
library(ggplot2)
str(tsne2List)

pt_size_sample <- 0.6; pt_alpha_sample <- 0.6
pt_size_cluster <- 0.6; pt_alpha_cluster <- 0.6
base_axis_size <- 20

pList_sample <- list(); pList_cluster <- list()
for(j in 1:10){ # length(tsne2List)
  message('j = ', j)
  pList_sample[[j]] <- plot_scatter(tsne2List[[j]], meta_data, label_name = 'sample', border_col='gray13',
                                    base_size=base_axis_size,palette_use=cols_sample, point_size = pt_size_sample, 
                                    point_alpha = pt_alpha_sample,no_guides=T) + labs(x=NULL, y=NULL)
  pList_cluster[[j]] <- plot_scatter(tsne2List[[j]], meta_data, label_name = 'cell_type', base_size=20,
                                     palette_use=cols_cluster2, point_size = pt_size_cluster, border_col='gray13',
                                     point_alpha =pt_alpha_cluster,no_guides=T) + labs(x=NULL, y=NULL)
}
names(pList_sample) <-names(pList_cluster) <-  names(tsne2List)
pList_sample[[1]]
library(cowplot)
subnames <- c("iDR-SC" ,"Seurat V3", "Harmony","fastMNN", "Uncorrected")

p12 <- plot_grid(plotlist =pList_sample[subnames], nrow=1, ncol=5)
ggsave(file="tsne_sample_heatmap_Hip2_merge70.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 1000)
pList_cluster[[1]]

p1 <- pList_cluster[[1]]
pList_cluster[[1]] <- Seurat::  LabelClusters(p1, id='cell_type', repel = T, 
                                              size = 6, color = 'white', 
                                              box = T, position = "nearest")
p12 <- plot_grid(plotlist =pList_cluster[subnames], nrow=1, ncol=5)
ggsave(file="tsne_cluster_heatmap_Hip2_merge70.png", plot = p12,
       width = 25, height =5, units = "in", dpi = 500)
plot_scatter(tsne2List[[1]], meta_data, label_name = 'cell_type', base_size=20,
             palette_use=cols_cluster2, 
             point_size = pt_size_cluster, 
             point_alpha =pt_alpha_cluster,no_guides=F)

# subnames2 <- c("Scanorama",  "scGen", "scVI", "MEFISTO", "PASTE")
# p12 <- plot_grid(plotlist =pList_sample[subnames2], nrow=1, ncol=5)
# ggsave(file="tsne_Py4_sample_heatmap_Hip2_merge70.png", plot = p12,
#        width = 25, height =5, units = "in", dpi = 200)
# p12 <- plot_grid(plotlist =pList_cluster[subnames2], nrow=1, ncol=5)
# ggsave(file="tsne_Py4_cluster_heatmap_Hip2_merge70.png", plot = p12,
#        width = 25, height =5, units = "in", dpi = 200)


##### 5c #####
rawID <- c(10, 11, 1, 7, 17, 8, 15, 13, 3, 9,  4,  16,2, 6, 5, 14, 12)
names(rawID) <- 1: 17
load("housekeep_dat_degs_Hip2_merge70.rds")
dat_degs$cluster <- factor(replace_ID(as.numeric(dat_degs$cluster), rawID = rawID), 
                           levels = 1:length(rawID))

dat_degs <- dat_degs[order(dat_degs$cluster),]
load("seuAll_hip2_merge70.rds")
library(dplyr)
library(Seurat)
n <- 5
dat_degs %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC) -> top10

seu@assays$RNA@var.features <- row.names(seu)
seu <- ScaleData(seu)
seus <- subset(seu, downsample = 4000)
color_id <- as.numeric(levels(Idents(seus)))


## HeatMap
p1 <- doHeatmap(seus, features = top10$gene, cell_label= "Domain",
                grp_label = F, grp_color = cols_cluster_low_reso[color_id],
                pt_size=6,slot = 'scale.data') + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=7, face= "italic", family='serif'))

ggsave(paste0('./Hip2All_merge70',"_top",n,"DEGs_heatmap_reOrder.pdf"), plot = p1, 
       width = 10, height = 8, units = "in", dpi = 1000)
##### 5d #####
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
load("Hip2_Merge70_clusterK17_renumber_idrsc.rds")
st_data <- readRDS("mergedSquare70_seulist_slideV2_mouseHip2.RDS")

## combine the weight matrix first
slice_set <- 1
norm_weights = read.csv(paste('./hippo_decon_update/output_weights_my_dge_hippo_mergedSquare70_st',slice_set,'.csv',sep=""),row.names=1)
norm_weights[1:4, 1:5]
dim(norm_weights)
row.names(norm_weights) <- paste0(row.names(norm_weights), "_hip1")
slice_set <- 2
norm_weights2 = read.csv(paste('./hippo_decon_update/output_weights_my_dge_hippo_mergedSquare70_st',slice_set,'.csv',sep=""),row.names=1)
norm_weights2[1:4, 1:5]
row.names(norm_weights2) <- paste0(row.names(norm_weights2), "_hip2")
norm_weights <- rbind(norm_weights, norm_weights2)
dim(norm_weights)
colnames(norm_weights)

## revise cell names
colnames(norm_weights)[colnames(norm_weights) =='Choroid_plexus'] = 'Choroid plexus'
colnames(norm_weights)[colnames(norm_weights) =='Cajal_Retzius'] = 'Cajal-Retzius'
colnames(norm_weights)[colnames(norm_weights) =='Endothelial_stalk'] = 'Endothelial stalk'
colnames(norm_weights)[colnames(norm_weights) =='Endothelial_tip'] = 'Endothelial tip'
colnames(norm_weights)[colnames(norm_weights) =='Entorhinal.cortex'] = 'Entorhinal cortex'
colnames(norm_weights)[colnames(norm_weights) =='Deep.layer.subiculum'] = 'Deep layer subiculum'
colnames(norm_weights)[colnames(norm_weights) =='Entorhinal.cortex..IEG.'] = 'Entorhinal cortex (IEG)'
colnames(norm_weights)[colnames(norm_weights) =='Dentate.principal'] = 'Dentate principal'
colnames(norm_weights)[colnames(norm_weights) =='CA1.principle.anterior'] = 'CA1 principal (Anterior)'
colnames(norm_weights)[colnames(norm_weights) =='CA1.principal'] = 'CA1 principal'
colnames(norm_weights)[colnames(norm_weights) =='Anterior.Subiculum..proximal.to.CA1'] = 'Anterior subiculum, proximal to CA1'

colnames(norm_weights)[colnames(norm_weights) =='Lateral.CA3.principal'] = 'CA3 principal (Lateral)'
colnames(norm_weights)[colnames(norm_weights) =='Dentate.hilum'] = 'Dentate hilum'
colnames(norm_weights)[colnames(norm_weights) =='MyelinProcesses'] = 'Myelin processes'
colnames(norm_weights)[colnames(norm_weights) =='Layer.6'] = 'Layer 6'
colnames(norm_weights)[colnames(norm_weights) =='Layer.5'] = 'Layer 5'
colnames(norm_weights)[colnames(norm_weights) =='Medial.habenula'] = 'Medial habenula'
colnames(norm_weights)[colnames(norm_weights) =='Medial.habenular..lateral.portion'] = 'Medial habenula (Lateral)'

colnames(norm_weights)[colnames(norm_weights) =='Lateral.habenula'] = 'Lateral habenula'
colnames(norm_weights)[colnames(norm_weights) =='MyelinProcesses'] = 'Myelin processes'
colnames(norm_weights)[colnames(norm_weights) =='Lateral.Geniculate.Nucleus...LGN.'] = 'Lateral geniculate nucleus'
colnames(norm_weights)[colnames(norm_weights) =='Lateral.dorsal.thalamic.nuclei..LD.'] = 'Lateral dorsal thalamic nuclei (LD)'
####
level_names <- c("CA1 principal (Anterior)", "CA1 principal", "Anterior subiculum, proximal to CA1", 
                 "CA2", "CA3", "CA3 principal (Lateral)", "Dentate principal",  "Dentate hilum" , 
                 "Interneuron",  "Neuron"  ,  "Neurogenesis"  ,                                               
                 "Macrophage"  ,                        "Ependymal"  ,
                 "Layer 6"      ,                       "Layer 5"  , 
                 "Microglia"  , "Astrocyte" ,     "Oligodendrocyte" ,
                 "Choroid plexus",                                        
                 "Cajal-Retzius", "Endothelial stalk",  "Endothelial tip"   ,              
                 "Entorhinal cortex",    "Entorhinal cortex (IEG)"  ,             
                 "Deep layer subiculum"   ,    "Mural" ,                    
                 "Myelin processes"   ,                 "Polydendrocyte" ,                    
                 "Medial habenula"  ,                   "Medial habenula (Lateral)",   
                 "Lateral habenula" ,                   "Lateral geniculate nucleus"  ,       
                 'Lateral dorsal thalamic nuclei (LD)'   )


cluster_id = as.data.frame(unlist(cluster_idrsc_renumber))
rownames(cluster_id) = c(paste0(colnames(st_data[[1]]), "_hip1"), 
                         paste0(colnames(st_data[[2]]), "_hip2"))
colnames(cluster_id) = "Cluster"
cluster_weight = merge(norm_weights,cluster_id,by = 0)
rownames(cluster_weight) = cluster_weight$Row.names
cluster_weight = cluster_weight[,-1]

percentage = as.data.frame(cluster_weight %>% group_by(Cluster) %>% summarise(across(everything(), sum)))
percentage = percentage[order(as.numeric(percentage$Cluster)),]
percentage$Cluster = as.numeric(percentage$Cluster)
percentage_long <- melt(as.data.frame(percentage),id.vars ='Cluster')

cols_group <- c(ggthemes::tableau_color_pal("Tableau 20")(20), 
                ggthemes::tableau_color_pal("Classic 20")(20)[c(5,7,11, 15)],
                tune_colors(ggthemes::tableau_color_pal("Classic 10 Light")(10)[-1], alpha = 0.4))
cols_group <- c(cols_group[1:2],"#aec7e866", "tan2", cols_group[3:33])

scale_flag <- c(TRUE,FALSE)
need_scale <- T
pList <- list()
for(ii in 1:2){
  
  if (need_scale == scale_flag[[ii]]){
    geom_bar_position = 'stack'
  }else{
    geom_bar_position = 'fill'
  }
  if (need_scale == scale_flag[[ii]]){
    percentage_scale = matrix(0,nrow(percentage),ncol(percentage))
    for (tt in 1:nrow(percentage)){
      percentage_scale[tt,] = as.numeric(cbind(percentage$Cluster[tt],percentage[tt,2:ncol(percentage)]/sum(rowSums(percentage[,2:ncol(percentage)]))))
    }
    colnames(percentage_scale) = colnames(percentage)
    #sum(percentage_scale[,2:ncol(percentage)])
    percentage_long <- melt(as.data.frame(percentage_scale),id.vars ='Cluster')
  }
  percentage_long$variable <- factor(percentage_long$variable, levels=level_names)
  
  if(ii==1){
    pList[[ii]] <- bar_plot_hip2(percentage_long,geom_bar_position=geom_bar_position, 
                            legend_position='none', fill_color = cols_group,
                            title_name = paste0(""))
  }else{
    pList[[ii]] <- bar_plot_hip2(percentage_long,geom_bar_position=geom_bar_position, 
                            legend_position='right', fill_color = cols_group,
                            title_name = paste0(""))
  }
  pList[[ii]] <- pList[[ii]] + theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 14))
  
  
}
figure_celltype_percentage <- ggarrange(plotlist = pList,ncol=2,nrow=1,align='hv',common.legend = T,
                                        legend = 'right')

ggsave(file='./result/combined_Low_percentage.png', plot =figure_celltype_percentage, 
       width = 19, height = 5, units = "in", bg = 'white', dpi = 1000,limitsize = F)
##### 5e #####
##### Start plot the pseudotime distribution on spatial coordinates
load("pseudotime_tmpList_hip2_merge70.rds")

pList_pseudo <- list()
for(r in 1:2){ ## each sample
  # r <- 1
  if(r ==1) pt_size = 0.25
  if(r == 2) pt_size = 0.1
  message("r = ", r)
  pseudotime_tmp <- pseudotime_tmpList[[r]]  [-idx_outList[[r]]]
  indx <- which(!is.na(pseudotime_tmp))
  pseudotime_tmp <- pseudotime_tmp[indx]
  pos_tmp <- coordinate_rotate(posList_raw[[r]][-idx_outList[[r]],][indx,], rotate_angles[r])
  
  dat1 <- data.frame(row=pos_tmp[,1], col=pos_tmp[,2], pseudotime=pseudotime_tmp)
  med <- quantile(dat1$pseudotime, 0.3)
  ptmp <- ggplot(dat1, aes(x=row, y=col, color=pseudotime)) + 
    geom_point(size=pt_size, alpha=1) +
    #scale_color_gradientn(colours = c("red", "green"))+ 
    scale_color_gradientn(colours = c( "#FFBB78", "#96B5E9", "#F54902"))+
    #scale_color_gradientn(colours = c("#96B5E9", "#F54902"))+
    #scale_color_gradient2(low = "#00008BCC",mid="#228B22CC", high = "#FFFF33", midpoint = med)+
    theme_hip2()+   scale_x_reverse() + theme(legend.position = 'none')
  
  pList_pseudo[[r]]<- ptmp
  
  
}
library(cowplot)
#pList_pseudo[[2]]
p12 <- plot_grid(plotlist =pList_pseudo, nrow=1, ncol=2, byrow=T)
ggsave(file=paste0("Hip2_merge70_psedotime_spatialheatmapV2.png"), plot = p12,
       width = 7, height =3, units = "in", dpi = 400)



load("pseudo.slingshot_hip2_merge70.rds")
load("hip2_merge70_tsne2_allMethods.rds")
tsne_idrsc <- tsne2List[[1]]

dat1 <- data.frame(PC1=tsne_idrsc[,1], PC2=tsne_idrsc[,2], pseudotime=pseudo.slingshot,
                   scaled_pseudotime = range01(pseudo.slingshot, na.rm=T))
library(ggplot2)

med <- quantile(dat1$scaled_pseudotime, 0.3, na.rm=T)
ggplot(dat1, aes(x=PC1, y=PC2, color=scaled_pseudotime)) + geom_point(alpha=1) +
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
  cowplot::theme_cowplot()+ xlab('tSNE 1') + ylab("tSNE 2")
