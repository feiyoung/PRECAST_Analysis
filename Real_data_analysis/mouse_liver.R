library(PRECAST)
library(Seurat)

dir_url <- "https://github.com/feiyoung/PRECAST_Analysis/raw/main/Real_data_analysis/data/seulist_mouseLiverST8.RDS"

seulist <- readRDS(url(dir_url))
yList <- lapply(seulist, function(x) x$cluster)
table(unlist(yList)) ## K= 6
save(yList, file='yList_mouseLiver8.rds')

library(Seurat)
genelist <- row.names(seulist[[1]])
getXList <- PRECAST:::getXList
datList <- getXList(seulist, genelist)

AdjList <- pbapply::pblapply(datList$posList, DR.SC::getAdj_auto, radius.upper= 10)

XList <- lapply(datList$XList, scale, scale=F)



# Integration analysis using PRECAST --------------------------------------


## 
hq <- 15
K_set <- 3:11
num_core <- 5
tic <- proc.time() # 
set.seed(1)
resList <- ICM.EM(XList,posList=NULL, AdjList=AdjList, q=hq, K=K_set, maxIter=30,
                 Sigma_equal =F, verbose=T, coreNum = num_core, coreNum_int = 1)
toc <- proc.time()
(time_used <- toc[3] - tic[3]) 
save(resList, time_used, file='resList_idrsc_chooseK_mouseLiverST8.rds')

reslist <- SelectModel(resList)


clusterList <- reslist$cluster
save(clusterList, file='idrsc_cluster7_mouseLiver8.rds')
cluster_metric <- function(hy, y, type='ARI'){
  
  require(mclust)
  require(aricode)
  switch(type, 
         ARI= adjustedRandIndex(hy, y),
         NMI = NMI(as.vector(hy), y))
}
combine_metric <- c(ARI=cluster_metric(unlist(reslist$cluster), unlist(yList)),
 NMI=cluster_metric(unlist(reslist$cluster), unlist(yList), type="NMI"))
combine_metric
sep_metric <- sapply(1:8, function(j) cluster_metric((reslist$cluster[[j]]), (yList[[j]])))
sep_metric


posList <- datList$posList
save(posList, file="posList_mouseLiver8.rds")
matlist2mat <- PRECAST:::matlist2mat
hZ_idrsc <- matlist2mat(reslist$hZ)

set.seed(1)
hZtsne_precast <- scater::calculateTSNE(t(hZ_idrsc))


# Downstream Analysis -----------------------------------------------------


# Combined DEG analysis ---------------------------------------------------


data(Mouse_HK_genes, package='PRECAST')
head(Mouse_HK_genes)
## find the housekeeping genes used in iDRSC model
houseKeep <- intersect(genelist, Mouse_HK_genes$Gene )

save(houseKeep, file="housekeep_genes_mouseLiver8.rds")

reslist <- SelectModel(resList, return_para_est=T)
Rf <- reslist$Rf
get_correct_exp <- PRECAST:::get_correct_exp
hX <- get_correct_exp(XList, Rf, houseKeep)
colnames(hX) <- firstup(colnames(hX))

library(Seurat)
library(Matrix)
count <- sparseMatrix(i=1,j=1, x=1, dims=dim(t(hX)))
row.names(count) <- colnames(hX)
colnames(count) <- row.names(hX)
seuAll <- CreateSeuratObject(counts = count)
length(unique(colnames(seuAll)))
row.names(hX) <- colnames(seuAll)
seuAll[["RNA"]]@data <- t(hX)
Idents(seuAll) <- factor(unlist(reslist$cluster))

## Get the DEGs including negative log fold changes for Vacano plots
dat_degs_all <- FindAllMarkers(seuAll)
save(dat_degs_all, file='housekeep_dat_degsAll_mouseLiver8.rds')

load('housekeep_dat_degsAll_mouseLiver8.rds')

filename <- 'mouseLiver8_DEGsList_jointAnalysis_p0001.xlsx'
cutoff <- 0.001#

for(k in 1:7){
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs_all,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Domain', k),row.names = 0,
                   append = T)
}


n <- 10
dat_degs_all %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC) -> top10
seu <- seuAll
#seu <- NormalizeData(seu)
seu@assays$RNA@var.features <- row.names(seu)
seu <- ScaleData(seu)
seu[['RNA']]@data[1:4,1:4]
cols_cluster <- cols_cluster2
color_id <- as.numeric(levels(Idents(seu)))

genes_use <- c("Cyp2e1", "Oat", "Cyp2c37", "Gulo", "Glul", "Slc1a2", 'Cyp2d9', "Gstm3", 
               "Malat1", "Cox1","Hamp2", "Cyp3a44", "Gsn", "Dpt","Vim", "Tagln",  
               "Cyp2f2", "Sds", "Hal", "Ctsc", "Aldh1b1", "Hsd17b13", "Spp1")
which(! genes_use %in% top10$gene)
p1 <- doHeatmap(seu, features = genes_use, cell_label= "Domain",
                grp_label = F, grp_color = cols_cluster[color_id],
                #disp.max = 2.1,
                pt_size=6,slot = 'scale.data') + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=12, face= "italic", family='serif'))
ggsave(paste0("Liver8All_usedGenesDEGs_heatmapV2.pdf"), plot = p1, width = 8, height = 6, units = "in", dpi = 1000)



# Domain DEG enrichment analysis ------------------------------------------


library(gprofiler2)

termList_cluster <- list()
for(k in 1: 7){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs_all,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  que1 <- toupper(dat_degs_sub3$gene)
  gostres <- gost(query = que1,
                  organism = "mmusculus", correction_method='fdr')
  termList_cluster[[k]] <- gostres
}
save(termList_cluster, file='profiler_termList_DomainDEGs_mouseLiver8.rds')

load('profiler_termList_DomainDEGs_mouseLiver8.rds')


source_set <- c("GO:BP", 'KEGG', "GO:CC", "GO:MF", "HPA",  "REAC") # 
## KEGG
source1 <- source_set
ss <- which(source_set %in% source1)
ntop = 5
cols <- c("steelblue3", "goldenrod", "brown3", "#f98866", "#ff420e", "green4" )

names(cols) <- source_set
pList_enrich <- list()
for(ii in 1:7){
  ## ii <- 1
  message("ii=", ii)
  gostres2 <- termList_cluster[[ii]]
  max(gostres2$result$term_size)
  dat_tmp <- subset(gostres2$result, term_size<500)
  table(dat_tmp$source)
  dat1 <- get_top_pathway1(dat_tmp, ntop=ntop, source_set = source_set)
  #dat1 <- rbind(dat1, subset(dat_tmp[,c("source", "term_name", "term_id","p_value")]))
  dat1$nlog10P <- -log10(dat1$p_value)
  dat1_sub <- subset(dat1[order(dat1$nlog10P),], source %in% source1 & nchar(term_name)< 60)
  pList_enrich[[ii]] <- barPlot_enrich(dat1_sub, source='source', 'term_name',
                                       'nlog10P', cols=cols[source_set[ss]],
                                       base_size = 25, bar_width = 0.9) + ylab("-log10(p-adj)")+
    xlab("Biological terms") + ggtitle(paste0("Domain", ii))
  
}


p12 <- patchwork::wrap_plots(pList_enrich, nrow=4, ncol=2)
ggsave(file=paste0("Domain1_7_DEG_enrich_barplot_mouseLiverV2.png"), plot = p12,
       width = 28, height =26, units = "in", dpi = 200)


# Explore the subtype of veins' function ----------------------------------
k_set1 <- 1:2
central_vein1_genes <- subset(dat_degs_all,  cluster==k_set1[1] & p_val_adj<cutoff & avg_log2FC>0.25)$gene
central_vein2_genes <- subset(dat_degs_all,  cluster==k_set1[2] & p_val_adj<cutoff & avg_log2FC>0.25)$gene
unique_Cvein1 <- setdiff(central_vein1_genes, central_vein2_genes)
unique_Cvein2 <- setdiff(central_vein2_genes, central_vein1_genes)
library(gprofiler2)
que1 <- toupper(unique_Cvein1)
gostres_Cvein1 <- gost(query = que1, organism = "mmusculus", correction_method='fdr')
term_df_Cvein1 <- subset(gostres_Cvein1$result, term_size<500)

source_set <- c("GO:BP", 'KEGG', "GO:CC", "GO:MF", "HPA",  "REAC") # 
dat1 <- get_top_pathway1(term_df_Cvein1, ntop=5, source_set = source_set)
#dat1 <- rbind(dat1, subset(dat_tmp[,c("source", "term_name", "term_id","p_value")]))
dat1$nlog10P <- -log10(dat1$p_value)
dat1_sub <- subset(dat1[order(dat1$nlog10P),], source %in% source_set & nchar(term_name)< 60)
barPlot_enrich(dat1_sub, source='source', 'term_name',
               'nlog10P', cols=cols,
               base_size = 25, bar_width = 0.9) + ylab("-log10(p-adj)")+
  xlab("Biological terms") + ggtitle(paste0("Central vein subtype 1"))

que1 <- toupper(unique_Cvein2)
gostres_Cvein2 <- gost(query = que1, organism = "mmusculus", correction_method='fdr')
term_df_Cvein2 <- subset(gostres_Cvein2$result, term_size<500)

dat1 <- get_top_pathway1(term_df_Cvein2, ntop=5, source_set = source_set)
dat1$nlog10P <- -log10(dat1$p_value)
dat1_sub <- subset(dat1[order(dat1$nlog10P),], source %in% source_set & nchar(term_name)< 60)
barPlot_enrich(dat1_sub, source='source', 'term_name',
               'nlog10P', cols=cols,
               base_size = 25, bar_width = 0.9) + ylab("-log10(p-adj)")+
  xlab("Biological terms") + ggtitle(paste0("Central vein subtype 2"))


k_set1 <- 6:7
portal_vein1_genes <- subset(dat_degs_all,  cluster==k_set1[1] & p_val_adj<cutoff & avg_log2FC>0.25)$gene
portal_vein2_genes <- subset(dat_degs_all,  cluster==k_set1[2] & p_val_adj<cutoff & avg_log2FC>0.25)$gene
unique_Pvein1 <- setdiff(portal_vein1_genes, portal_vein2_genes)
unique_Pvein2 <- setdiff(portal_vein2_genes, portal_vein1_genes)

que1 <- toupper(unique_Pvein1)
gostres_Pvein1 <- gost(query = que1, organism = "mmusculus", correction_method='fdr')
term_df_Pvein1 <- subset(gostres_Pvein1$result, term_size<500)

dat1 <- get_top_pathway1(term_df_Pvein1, ntop=5, source_set = source_set)
dat1$nlog10P <- -log10(dat1$p_value)
dat1_sub <- subset(dat1[order(dat1$nlog10P),], source %in% source_set & nchar(term_name)< 60)
barPlot_enrich(dat1_sub, source='source', 'term_name',
               'nlog10P', cols=cols,
               base_size = 25, bar_width = 0.9) + ylab("-log10(p-adj)")+
  xlab("Biological terms") + ggtitle(paste0("portal vein subtype 1"))

que1 <- toupper(unique_Pvein2)
gostres_Pvein2 <- gost(query = que1, organism = "mmusculus", correction_method='fdr')
term_df_Pvein2 <- subset(gostres_Pvein2$result, term_size<500)
dat1 <- get_top_pathway1(term_df_Pvein2, ntop=5, source_set = source_set)
dat1$nlog10P <- -log10(dat1$p_value)
dat1_sub <- subset(dat1[order(dat1$nlog10P),], source %in% source_set & nchar(term_name)< 60)
barPlot_enrich(dat1_sub, source='source', 'term_name',
               'nlog10P', cols=cols,
               base_size = 25, bar_width = 0.9) + ylab("-log10(p-adj)")+
  xlab("Biological terms") + ggtitle(paste0("portal vein subtype 2"))


# Trajectory inference ----------------------------------------------------
domain_dup_names <- c("CV-1", "CV-2", "Endo", "Hep", "Mes", "PV-1", "PV-2")

domain_ID <- 1:7
names(domain_ID) <- domain_dup_names

replace_ID <- function(x, rawID){
  # x <- cluster_idrsc[[1]]
  uni_x <- unique(x)
  if(!all(uni_x %in% rawID) ) stop("Some element of x is not in rawID!")
  
  idx <- which(rawID %in% uni_x)
  sub_rawID <- rawID[idx]
  y <- rep(NA, length(x))
  
  for(k in 1:length(sub_rawID)){
    y[x==sub_rawID[k]] <- names(sub_rawID)[k]
  }
  return(y)
}


cl_idrsc_cellRegion <- replace_ID(unlist(clusterList), domain_ID)
table(cl_idrsc_cellRegion)
cl_idrsc_cellRegion2 <- cl_idrsc_cellRegion
cl_idrsc_cellRegion2[cl_idrsc_cellRegion=="CV-1" | cl_idrsc_cellRegion=="CV-2"] <- "CV"
cl_idrsc_cellRegion2[cl_idrsc_cellRegion=="PV-1" | cl_idrsc_cellRegion=="PV-2"] <- "PV"
library(slingshot)
rd <- hZList[[1]]
# Infer pseudotimes
sds <- slingshot(rd, cl_idrsc_cellRegion2) #



## Get pseudotimes
ptall <- slingPseudotime(sds)


pseudo.slingshot2 <-  rowMeans(slingPseudotime(sds), na.rm=TRUE)
aggregate(pseudo.slingshot2, by=list(cl_idrsc_cellRegion), mean)
save(pseudo.slingshot2, cl_idrsc_cellRegion, file='pseudotime_allspots_mouseLiver8.rds')



library(SingleCellExperiment)
sce.liver <- SingleCellExperiment(assays=list(logcounts= t(hX)))
sce.liver$domain <- factor(unlist(clusterList), levels=1:7)
save(sce.liver, file='Liver8_sce_changedGenes_along_pseudotime.rds')



library(TSCAN)
pseudo <- testPseudotime(sce.liver, pseudotime=pt)

save(pseudo, pt,file='Liver8_pseudo_TSCANtest.rds')

#### Plot trajectory heatmap###############
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
features_reorder <- c("Glul", "Oat","Slc1a2","Cyp2e1","Cyp2a5",
                      "Gulo", "Cyp2c37","Lect2" , "Cyp2c29","Cyp2c69", "Cyp3a11","Rnase4", "Aldh1a1",
                      ## lower
                      "Cyp2f2", "Hal" ,"Sds", "Hsd17b13", "Ctsc","Rida", "A1bg")
#"Mup3",  "Hpx", "Scd1","Car3")

p1 <- scater::plotHeatmap(sce.liver, order_columns_by="Pseudotime", 
                          colour_columns_by= "domain", features= features_heat,
                          center=TRUE, swap_rownames="SYMBOL", cluster_rows=F)
dir.file <- ''
ggsave(paste0(dir.file, 'hcc4_Joint_Traject_HeatMap.pdf'), plot = p1, width = 12, height = 5, units = "in", dpi = 1000)




# Cell-type deconvolution -------------------------------------------------
## Reference data: http://bis.zju.edu.cn/MCA/gallery.html?tissue=Adult-Liver

library(spacexr)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
# library(ggpubr)
library(ggthemes)
library(cowplot)
library(Matrix)
library(data.table)
library(reshape2)

######## IMPORT DATA
###ST DATA

RCTD_structure <- function(sc_obj, clust_vr) {
  
  sc_obj[["Name"]] = sc_obj@meta.data[, clust_vr]
  
  # Cell type dictionary between cluster and cell-type
  ct <- unique(sc_obj@meta.data[, clust_vr])
  df_ct <- data.frame("Cluster" = 1:length(ct),
                      "Name" = ct)
  metadata <- sc_obj@meta.data %>%
    # Rownames to columns must be before left join since after it the rownames are erased
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(df_ct, by = c("Name" = "Name")) %>%
    # Change names to 鈥渂arcode鈥? 鈥渃luster鈥? 鈥渘UMI鈥?   
    mutate(
      cluster = Cluster,
      nUMI = nCount_RNA
    ) %>%
    dplyr::select(barcode, cluster, nUMI)
  
  expr_mtrx <- sc_obj@assays$RNA@counts
  
  return(list("meta_data" = metadata,
              "cell_type_dict" = df_ct,
              "dge" = expr_mtrx))
}


###SC DATA

meta_dat <- data.table::fread(file='./Deconvolution_dat/Adult-Liver_barcodes_anno.csv')

exp_data <- data.table::fread(file='./Deconvolution_dat/Adult-Liver_dge.csv')
gene_data <- data.table::fread(file='./Deconvolution_dat/Adult-Liver_gene.csv', header = T)
counts <- Seurat::as.sparse(exp_data[,-1])
row.names(counts) <- gene_data$x

meta_dat <- as.data.frame(meta_dat)
row.names(meta_dat) <- meta_dat$V1
setdiff(row.names(meta_dat), colnames(counts))
setdiff(colnames(counts), row.names(meta_dat))
sc_data <- Seurat::CreateSeuratObject(counts=counts, meta.data= meta_dat)
saveRDS(sc_data, file='./Deconvolution_dat/sc_data_mouseLiver8_reference.RDS')
sc_data <- readRDS(file='./Deconvolution_dat/sc_data_mouseLiver8_reference.RDS')
colnames(sc_data@meta.data)[5] <- "celltype"
cell_type_columns <- "celltype"
table(sc_data@meta.data[,cell_type_columns])
delete_cell_type <- c("Epithelial cell_Lcn2 high", "Neutrophil_Ngp high", "Hepatocyte_mt high")
cell_used <- !(sc_data@meta.data[,cell_type_columns] %in% delete_cell_type)
sc_data_sub <-sc_data[, cell_used]
table(sc_data_sub@meta.data[,cell_type_columns])
for (i in 1: 8){
  # i <- 1
  message("i = ", i)
  hcc_st_s = seulist[[i]]
  st_count = hcc_st_s@assays$RNA@counts
  st_coord = hcc_st_s@meta.data[,c('row','col')]
  output_name = paste('mouseLiver8_st',i,sep="")
  output_dir = 'result_new'
  RCTD_run(sc_data_sub,cell_type_columns,st_count,st_coord,output_dir,output_name)   
}