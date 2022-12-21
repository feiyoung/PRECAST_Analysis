### dir_current <- "F:/Research paper/IntegrateDRcluster/AnalysisCode/PRECAST_Analysis/"
dir_current <- "./"
source(paste0(dir_current, "Real_data_results/Rcode/utility_funcs.R"))
##### 2a #####
setwd(paste0(dir_current, "Real_data_results/dataFiles/brain12/") )

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


##### 2h #####
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


