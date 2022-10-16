library(PRECAST)
library(Seurat)
segment_square <- function(pos, sq_nspots=70, by_order=T, verbose=T){
  tmp <- pos[,1]
  if(by_order){
    x_cut <- sort(tmp)[seq(1, length(tmp), length.out=sq_nspots+1)]
  }else{
    x_cut <- seq(min(tmp), max(tmp), length.out=sq_nspots+1)
  }
  
  tmp <- pos[,2]
  if(by_order){
    y_cut <- sort(tmp)[seq(1, length(tmp), length.out=sq_nspots+1)]
  }else{
    y_cut <- seq(min(tmp), max(tmp), length.out=sq_nspots+1)
  }
  
  i <- 1
  pos_new <- matrix(NA, sq_nspots^2, 2)
  areaList <- list()
  for(i1 in 1:sq_nspots){
    if(verbose)
      message('i1 = ', i1)
    for(i2 in 1:sq_nspots){
      if(i1 < sq_nspots && i2 < sq_nspots){
        tmp <- which(x_cut[i1] <=pos[,1] & pos[,1]<x_cut[i1+1] & y_cut[i2] <= pos[,2]& pos[,2] < y_cut[i2+1])
      }else if(i1 < sq_nspots && i2 == sq_nspots){
        tmp <- which(x_cut[i1] <=pos[,1] & pos[,1]< x_cut[i1+1] & y_cut[i2] <= pos[,2]& pos[,2] <= y_cut[i2+1])
      }else{
        tmp <- which(x_cut[i1] <=pos[,1] & pos[,1]<=x_cut[i1+1] & y_cut[i2] <= pos[,2]& pos[,2] < y_cut[i2+1])
      }
      
      areaList[[i]] <- tmp
      pos_new[i, ] <- c((x_cut[i1]+ x_cut[i1+1])/2, (y_cut[i2]+y_cut[i2+1])/2 )
      i <- i + 1
    }
  }
  idx <- which(sapply(areaList, function(x) length(x)>0))
  return(list(spotID_list=areaList[idx], pos_new = pos_new[idx, ]))
}

get_merged_seu <- function(seu, areaList, pos_new){
  require(Seurat)
  n_area <- length(areaList)
  count_new <- matrix(NA, nrow(seu), n_area)
  colnames(count_new) <- paste0("merge_spot", 1:n_area)
  row.names(count_new) <- row.names(seu)
  DefaultAssay(seu) <- "RNA"
  for(i in 1:n_area){ # 
    message('i = ', i, '/', n_area)
    if(length(areaList[[i]])>1){
      count_new[, i] <- rowSums(seu[["RNA"]]@counts[,areaList[[i]]])
    }else{
      count_new[, i] <- seu[["RNA"]]@counts[,areaList[[i]]]
    }
    
  }
  rm(seu)
  meta_data <- data.frame(row=pos_new[,1], col=pos_new[,2])
  row.names(meta_data) <- colnames(count_new)
  CreateSeuratObject(counts= as.sparse(count_new), meta.data = meta_data) 
}


barcodeList <- pbapply::pblapply(seuList_raw, function(x) colnames(x))
save(barcodeList, file='barcodeList_before_merge70_Bul20rep2.rds')
posList_before_merge70 <- pbapply::pblapply(seuList_raw, function(x) cbind(x$row, x$col))



#Segement to 70*70 #############
IDmap70List <- list()
seuList_square70 <- list()
for(r in 1:20){
  # r <- 1
  message("r = ", r)
  res_pos_seg1 <- segment_square(posList_before_merge70[[r]])
  IDmap70List[[r]] <- res_pos_seg1$spotID_list
  seu_new1 <- get_merged_seu(seuList_raw[[r]], res_pos_seg1$spotID_list, res_pos_seg1$pos_new)
  seuList_square70[[r]] <- seu_new1
}

saveRDS(seuList_square70, file='mergedSquare70_seuList_Bulb20_rep2.RDS')

sample_high_quality <- 1:16
barcodeList16 <- barcodeList[sample_high_quality]
save(barcodeList16, file='barcodeList16_before_merge70_Bulb16rep2.rds')
posList16_before_merge70 <- posList_before_merge70[sample_high_quality]
save(posList16_before_merge70, file='posList16_before_merge70_Bulb16rep2.rds')


IDmap70List16 <- IDmap70List[sample_high_quality]

length(unlist(res_pos_seg1$spotID_list))
save(IDmap70List16, file='IDmap70List16_Bulb16_rep2.rds')



seuList_square70 <- readRDS(file='mergedSquare70_seuList_Bulb20.RDS')
sample_high_quality <- 1:16
seuList16_merge70 <- seuList_square70[sample_high_quality]
## select top 2000 SVGs
seuList <- lapply(seuList16_merge70, DR.SC::FindSVGs, nfeatures=2000,num_core=1, verbose=TRUE)
spaFeatureList <- lapply(seuList, DR.SC::topSVGs, ntop=2000)

selectIntFeatures <- PRECAST:::selectIntFeatures
genelist <- selectIntFeatures(seuList,spaFeatureList)
filter_spot <- function(seu, min_feature=0){ # each spots at least include 10 non-zero features
  subset(seu, subset = nFeature_RNA > min_feature)
}

seulist <- pbapply::pblapply(seuList, function(x) x[genelist, ])
index_zeroList <- pbapply::pblapply(seuList, function(x) which(x$nFeature_RNA <= 15)) # empty set
seulist <- pbapply::pblapply(seulist, filter_spot)

saveRDS(seulist, file='slideV2_MOB16_rep2_merge70_seulist.RDS')
save(genelist, file='genelist_MOB16_rep2_merge70.rds')

# Integration analysis using PRECAST --------------------------------------

posList <- lapply(seulist, function(x) cbind(x$row, x$col))
genelist <- row.names(seulist[[1]])
getXList <- PRECAST:::getXList
datList <- getXList(seulist, genelist)


AdjList <- pbapply::pblapply(datList$posList, DR.SC::getAdj_auto, 
                             lower.med=4, upper.med=6, radius.upper= 90)


hq <- 15
K_set <- 2:14
num_core <- 10
tic <- proc.time() # 837058 spots
set.seed(1)
resList_merge70 <- ICM.EM(datList$XList,posList=NULL, AdjList=AdjList, 
                         q=hq, K=K_set, int.model=NULL,  beta_grid= seq(1,6, by=0.2), maxIter=30,
                         Sigma_equal =F, verbose=T, coreNum = num_core, coreNum_int = num_core)
toc <- proc.time()
(time_used <- toc[3] - tic[3]) 
lapply(resList_merge70[[1]]$cluster, table)
save(resList_merge70, time_used, file='resList_iDRSC_merge70_Bulb16_rep2.rds')


load("resList_iDRSC_merge70_Bulb16_rep2.rds")
selectModel(resList_merge70, pen_const = 1)$icMat
reslist <- selectModel(resList_merge70,  return_para_est = T)
attr(reslist, 'fit')$beta
clusterList <- reslist$cluster

save(clusterList, file = 'iDRSC_cluster12_Bulb18_merge70_rep2.rds')


## Downstream Analysis -----------------------------------------------------
### Combined DEG analysis   #####################################

Rf <- attr(reslist, "fit")$Rf
hX <- get_correct_exp(datList$XList, Rf, houseKeep)

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

save(seuAll, file='iDRSCbatchCorrected_seuAll_Bulb16_rep2_merge70.rds')

dat_degs <- FindAllMarkers(seuAll)
save(dat_degs, file='./housekeep_dat_degs_Bulb16_rep2_merge70.rds')

library(dplyr)
n <- 10
dat_degs %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC) -> top10
seu <- seuAll
#seu <- NormalizeData(seu)
seu@assays$RNA@var.features <- row.names(seu)
seu <- ScaleData(seu)
seu[['RNA']]@data[1:4,1:4]
seus <- subset(seu, downsample = 4000)
color_id <- as.numeric(levels(Idents(seus)))
seus[['RNA']]@data[1:4,1:4]
sum(top10$gene %in%  row.names(seu))
cols_cluster
cols_cluster2 <- cols_cluster
cols_cluster2[4] <- "green2"
## HeatMap
p1 <- doHeatmap(seus, features = top10$gene, cell_label= "Domain",
                grp_label = F, grp_color = cols_cluster2[color_id],
                pt_size=6,slot = 'scale.data') + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=4, face= "italic"))
ggsave(paste0('./Bulb16All_rep2_merge70',"_top",n,"DEGs_heatmap_Reorder.pdf"), plot = p1, 
       width = 10, height = 8, units = "in", dpi = 1000)

load("housekeep_dat_degs_Bulb16_rep2_merge70.rds")
dat_degs$cluster <- factor(replace_ID(as.numeric(dat_degs$cluster), rawID = rawID), levels = 1:length(rawID))


cutoff <- 0.001
dat_degs_sub <- subset(dat_degs,  avg_log2FC >0.25 &   p_val_adj<cutoff)
table(dat_degs_sub$cluster)
sum(table(dat_degs_sub$cluster))


filename <- 'Bulb16_rep2_merge70_DEGsList_jointAnalysis_Reordered001.xlsx'

for(k in 1: 12){
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Domain', k),
                   append = T)
}


library(gprofiler2)

termList_cluster <- list()
for(k in 1: 12){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  que1 <- toupper(dat_degs_sub3$gene)
  gostres <- gost(query = que1,
                  organism = "mmusculus", correction_method='fdr')
  termList_cluster[[k]] <- gostres
}
save(termList_cluster, file='profiler_termList_DomainDEGs_Bulb16_merge70.rds')

load('profiler_termList_DomainDEGs_Bulb16_merge70.rds')


source_set <- c("GO:BP","GO:CC", "GO:MF",   'KEGG')
## KEGG
source1 <- c("GO:BP","GO:CC", "GO:MF", "KEGG")
ss <- which(source_set %in% source1)
ntop = 4
cols <- c("steelblue3", "goldenrod", "brown3", "#f98866", "#ff420e" )
names(cols) <- source_set
pList_enrich <- list()
for(ii in 1:12){
  ## ii <- 1
  message("ii=", ii)
  gostres2 <- termList_cluster[[ii]]
  max(gostres2$result$term_size)
  dat_tmp <- subset(gostres2$result, term_size<500)
  table(gostres2$result$source)
  dat1 <- get_top_pathway1(dat_tmp, ntop=ntop, source_set = source_set)
  dat1$nlog10P <- -log10(dat1$p_value)
  
  pList_enrich[[ii]] <- barPlot_enrich(subset(dat1[order(dat1$nlog10P),], source %in% source1), source='source', 'term_name',
                                       'nlog10P', cols=cols[source_set[ss]],
                                       base_size = 21, bar_width = 0.9) + ylab("-log10(p-adj)")+
    xlab("Biological terms") + ggtitle(paste0("Domain", ii))
  
}

p12 <- patchwork::wrap_plots(plotlist = pList_enrich[1:12], nrow=6, ncol=2)

### Deconvolution analysis########################################


dat_raw <- data.table::fread(file="GSE121891_OB_6_runs.raw.dge.csv")

library(Matrix)
getSparseExpMat <- function(exprs){
  expMat <- exprs[, -1]
  expMat <- Seurat::as.sparse(expMat)
  #colnames(expMat) <- exprs[1,-1]
  row.names(expMat) <- as.vector(as.data.frame(exprs[,1])[,1]) # how to access to data.table.
  return(expMat)
}
expMat <- getSparseExpMat(dat_raw)
expMat[1:4,1:4]


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

source("AppRCTD_05_22.R")

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

meta_dat_fi2 <- data.table::fread(file='GSE121891_Figure_2_metadata.txt')

meta_data <- data.table::fread(file='GSE121891_OB_metaData_seurat.csv')

meta_data <- as.data.frame(meta_data)
row.names(meta_data) <- meta_data[,1]
sum(colnames(expMat) != row.names(meta_data))
expMat_inMeta <- expMat[,row.names(meta_data)]
barcode_all <- colnames(expMat_inMeta)
getsample_type <- function(x){
  sapply(x, function(y) strsplit(y, split = '_')[[1]][1])
}
sample_types <- getsample_type(barcode_all)
table(sample_types)
#  OC1   OC2   TR1   TR2   WT1   WT2
# 7288  9154  6705 10570  9997  7712
sample_types_fig2 <- getsample_type(meta_dat_fi2$V1)
table(sample_types_fig2)
### select the wild type data
indx_WT <- which(sample_types %in% c("WT1", "WT2"))

# expMat_WT <- expMat_inMeta[,indx_WT]
meta_WT <- meta_data[indx_WT, ]
table(meta_WT$ClusterName)
## Determine the neural subtype in the WT data meta_WT
neual_abbr <- paste0("N", 1:16)
index_neural <- which(meta_WT$ClusterName %in% neual_abbr)
neural_large_type <- meta_WT[index_neural, ]
index_WT_fig2 <- which(sample_types_fig2 %in% c("WT1", "WT2"))
meta_WT_fig2 <- as.data.frame(meta_dat_fi2[index_WT_fig2, ])
row.names(meta_WT_fig2) <- meta_WT_fig2$V1
table(meta_WT_fig2$FinalIds)
meta_WT_fig2$ClusterName_copy <- meta_WT_fig2$ClusterName
meta_WT_fig2$ClusterName <- meta_WT_fig2$FinalIds

col_used_names <- c('orig.ident', 'experiment', 'ClusterName')
data_meta_WT <- rbind(meta_WT[-index_neural, col_used_names],
                      meta_WT_fig2[,col_used_names])



sc_data = Seurat::CreateSeuratObject(counts = expMat_inMeta[, row.names(data_meta_WT)],
                                     meta.data = data_meta_WT)
head(sc_data)
saveRDS(sc_data, file='sc_data_MOB_reference_new.RDS')

cell_type_columns <- "ClusterName"

seulist <- readRDS(file="slideV2_MOB16_rep2_merge70_seulist.RDS")

for (i in 2: 16){
  # i <- 1
  message("i = ", i)
  hcc_st_s = seulist[[i]]
  st_count = hcc_st_s@assays$RNA@counts
  st_coord = hcc_st_s@meta.data[,c('row','col')]
  output_name = paste('bulb_st',i,sep="")
  output_dir = 'result_new'
  RCTD_run(sc_data,cell_type_columns,st_count,st_coord,output_dir,output_name)   
}

## Trajectory inference ################################
### Trajectory inference
library(SingleCellExperiment)
# Load batch corrected expression matrix
load("iDRSCbatchCorrected_seuAll_Bulb16_rep2_merge70.rds")
sce_traject <- SingleCellExperiment(assays=list(logcounts=seuAll[["RNA"]]@data))
### Low-dimensional embeddings from iDR-SC
load("iDRSC_lowdim_embedings_bulb16_rep2_merge70.rds")
rd <- hZ_idrsc 

### Domain clusters from iDR-SC
load("Bulb16_rep2_Merge70_clusterK12_renumber_idrsc.rds")
cl_drsc <- unlist(cluster_idrsc_renumber)
### transfer the domain IDs to cell type regions
group_number <- c(1, 1, 1, 1, 1, 2, 1, 4)
regions <- c("rostral migratory stream", "granule cell layer", "GCL/inner plexiform layer", 
             "mitra layer", "outer plexiform layer", "glomerular layer", "olfactory nerve layer", 
             "Low quality region")

domain_dup_names <- NULL
for(kk in 1: length(group_number)){
  
  if(group_number[kk]>1){
    tmp <- rep(regions[kk], group_number[kk])
  }else{
    tmp <- regions[kk]
  }
  domain_dup_names <- c(domain_dup_names, tmp)
}

domain_ID <- 1:12
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


cl_idrsc_cellRegion <- replace_ID(cl_drsc, domain_ID)
sce_traject$cellRegion <- cl_idrsc_cellRegion
reducedDim(sce_traject, "PCA") <- rd
## select the domains for trajectory inference

idx <- which(!(cl_idrsc_cellRegion %in% "Low quality region") )
sce_traject_sub <- sce_traject[,idx]
rd_sub <- rd[idx, ]
cl_idrsc_cellRegion_sub <- cl_idrsc_cellRegion[idx]
sce_traject_sub$cellRegion <- cl_idrsc_cellRegion_sub
set.seed(1)
library(slingshot)
# Infer pseudotimes
sds <- slingshot(rd_sub, cl_idrsc_cellRegion_sub) #
ptall <- slingPseudotime(sds)

## Check the pseudotime
pseudo.slingshot <-  rowMeans(slingPseudotime(sds), na.rm=TRUE)

pseudotime_allspots <- rep(NA, length(cl_idrsc_cellRegion))
pseudotime_allspots[idx] <- pseudo.slingshot

save(pseudotime_allspots, file='pseudotime_allspots_bulb16_merge.rds')