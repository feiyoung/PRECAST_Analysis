

library(DR.SC)
library(Seurat)
library(Matrix)
library(dplyr)
library(PRECAST)

url_hcc <- "https://github.com/feiyoung/PRECAST_Analysis/raw/main/Real_data_analysis/data/"
posList <- list()
seuList <- list()
for(iter in 1: 4){
  # iter <- 1
  message("iter = ", iter)
  hcc <- readRDS(paste0(url_hcc,"HCC", iter, "_seu.RDS"))
  seuList[[iter]] <- hcc
  posList[[iter]] <- cbind(row=hcc$row, col=hcc$col)
}

### Save seuList for HCC4 data##################
saveRDS(seuList, file='HCC4_seuList.RDS')
save(posList, file = 'HCC4_posList.rds')

load('HCC4_spatialFeatureList.rds')

selectIntFeatures <- PRECAST:::selectIntFeatures
getXList <- PRECAST:::getXList

genelist <- selectIntFeatures(seuList, spaFeatureList)
length(unique(genelist))
save(genelist, file='HCC4_genelist2000.rds')


datList <- getXList(seuList, genelist)
saveRDS(datList, file='HCC4_datList.RDS')

metadataList <- lapply(seuList, function(x) x@meta.data)
save(metadataList, file='metadataList_hcc4.rds')


XList <- datList$XList
posList <- datList$posList
indexList <- datList$indxList

## Integration analysis using PRECAST ############################################################
q <- 15; K <- 9# 2:11
# Integrate 12 samples
tic <- proc.time() # 
set.seed(1)
resList <- ICM.EM(XList, posList=posList, q=q, K=K, 
                 platform = 'Visium',maxIter = 30,
                 Sigma_equal =F, verbose=T, coreNum = 5, coreNum_int = 5)
toc <- proc.time()
time_used_chooseK <- toc[3] - tic[3] 

save(time_used_chooseK, resList, file ="idrsc_HCC4_chooseK.rds")  
reslist <- SelectModel(resList)

## Export this function
matlist2mat <- PRECAST:::matlist2mat

hZ_idrsc_K9 <- matlist2mat(reslist$hZ)
set.seed(1)
hZ_tsne_idrsc_K9 <- scater::calculateTSNE(t(hZ_idrsc_K9))

# Downstream analysis -----------------------------------------------------


## Combined DEG analysis#########################

load("Housekeeping_GenesHuman.RData")
head(Housekeeping_Genes)
## find the housekeeping genes used in iDRSC model
houseKeep <- intersect(genelist, Housekeeping_Genes$Gene.name )
save(houseKeep, file="housekeep_genes_HCC4.rds")

reslist <- selectModel(resList, return_para_est=T)
Rf <- attr(reslist, "fit")$Rf
hX <- get_correct_exp(XList, Rf, houseKeep)

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

Idents(seuAll) <- factor(unlist(cluster_idrsc_renumber))

dat_degs <- FindAllMarkers(seuAll, only.pos = T)
save(dat_degs, file='housekeep_dat_degs_HCC4.rds')

n <- 10
dat_degs %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC) -> top10
seu <- seuAll
#seu <- NormalizeData(seu)
seu@assays$RNA@var.features <- row.names(seu)
seu <- ScaleData(seu)
seu[['RNA']]@scale.data
seu[['RNA']]@data[1:4,1:4]
seus <- subset(seu, downsample = 4000)
color_id <- as.numeric(levels(Idents(seus)))

cols_cluster <- c("#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1",
                  "#1f77b4", "#A65628",  "blue", "green2")

## First letter use upper capital
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

counts <- seus[['RNA']]@counts
row.names(counts) <- firstup(row.names(count))
seus2 <- CreateSeuratObject(counts = counts)
data1 <- seus[['RNA']]@scale.data
row.names(data1) <- firstup(row.names(data1))
seus2[['RNA']]@scale.data <- data1
Idents(seus2) <- Idents(seus)
p1 <- doHeatmap(seus2, features = firstup(top10$gene), cell_label= "Domain",
                grp_label = F, grp_color = cols_cluster[color_id],
                pt_size=6,slot = 'scale.data') + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=7, face= "italic", family='serif'))
ggsave(paste0('HCCAll',"_top",n,"DEGs_heatmap_captalfirstletter.pdf"), plot = p1, width = 10, height = 8, units = "in", dpi = 1000)


row.names(dat_degs) <- firstup(row.names(dat_degs))
dat_degs$gene <-  firstup(dat_degs$gene)
filename <- 'HCC4_DEGsList_jointAnalysis_Reordered001.xlsx'
cutoff <- 0.001

for(k in 1:9){
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  
  dat <- as.data.frame(dat_degs_sub3)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Domain', k),
                   append = T)
}


### Domain DEG enrichment analysis#########################

library(gprofiler2)

termList_cluster <- list()
for(k in 1: 9){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs,  cluster==k & p_val_adj<cutoff & avg_log2FC>0.25)
  que1 <- toupper(dat_degs_sub3$gene)
  gostres <- gost(query = que1,
                  organism = "hsapiens", correction_method='fdr')
  termList_cluster[[k]] <- gostres
}
save(termList_cluster, file='profiler_termList_DomainDEGs_HCC4.rds')
source_set <- c("GO:BP", 'KEGG',  "REAC") # "GO:CC", "GO:MF", "HPA",
## KEGG
source1 <- source_set
ss <- which(source_set %in% source1)
ntop = 4
cols <- c("steelblue3", "goldenrod", "brown3", "#f98866", "#ff420e" )
extra_term_idList <- list()
extra_term_idList[[1]] <- c("KEGG:05207")
for(ii in 2:9) extra_term_idList[[9]] <- "NA"
extra_term_idList[[5]] <- c("KEGG:05207")
extra_term_idList[[4]] <- c("REAC:R-HSA-9656223", 'REAC:R-HSA-6802949')
names(cols) <- source_set
pList_enrich <- list()
for(ii in 1:9){
  ## ii <- 5
  message("ii=", ii)
  gostres2 <- termList_cluster[[ii]]
  max(gostres2$result$term_size)
  dat_tmp <- subset(gostres2$result, term_size<500)
  table(dat_tmp$source)
  dat1 <- get_top_pathway1(dat_tmp, ntop=ntop, source_set = source_set)
  dat1 <- rbind(dat1, subset(dat_tmp[,c("source", "term_name", "term_id","p_value")], term_id %in% extra_term_idList[[ii]]))
  dat1$nlog10P <- -log10(dat1$p_value)
  dat1_sub <- subset(dat1[order(dat1$nlog10P),], source %in% source1 & nchar(term_name)< 60)
  pList_enrich[[ii]] <- barPlot_enrich(dat1_sub, source='source', 'term_name',
                                       'nlog10P', cols=cols[source_set[ss]],
                                       base_size = 25, bar_width = 0.9) + ylab("-log10(p-adj)")+
    xlab("Biological terms") + ggtitle(paste0("Domain", ii))
  
}


p12 <- patchwork::wrap_plots(plotlist = pList_enrich, nrow=5, ncol=2)

## Correlation for different domains ---------------------------------------

hZ_idrsc <- hZList[[1]]
cluster_idrsc_vec <- as.numeric(unlist(cluster_idrsc_renumber))
subsample <- function(hZ_idrsc, cluster_idrsc_vec, sample_rate = 0.5){
  nn <- length(cluster_idrsc_vec)
  idx <- sort(sample(nn, nn* sample_rate))
  return(list(hZ=hZ_idrsc[idx, ], cluster=cluster_idrsc_vec[idx]))
}
sublist <- subsample(hZ_idrsc, cluster_idrsc_vec, 0.15)
order_idx <- order(sublist[[2]])
cluster_orderd <- sublist[[2]][order_idx]
corMat <- cor(t(sublist[[1]][order_idx, ]))
row.names(corMat) <- paste0('sp', 1:nrow(corMat))
colnames(corMat) <- paste0('spot', 1:nrow(corMat))


p1 <- doHeatmap.matrix(corMat, cluster_orderd, legend_title = 'Domain', grp_color=cols_cluster)



## SVGs test for each data ##########################
hZ_idrsc <- hZList[[1]]
SVGsList <- list()
for(r in 1:4){
  res_tmp <- conditionalSVGs(seulist[[r]], hZ_idrsc[indexList[[r]], ], genelist=genelist, num_core = 1)
  SVGsList[[r]] <- res_tmp
  save(SVGsList, file='HCC4_SVGs_SVGsList.rds') # use this method to savd data in time
}



# Deconvolution analysis --------------------------------------------------

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

######## IMPORT DATA
###ST DATA
hcc_st_list = readRDS("HCC4_seuList.RDS")

###SC DATA
meta_ori = read.csv("gse125449/GSE125449_Set1_samples.txt",sep='\t')
sc_count = as.matrix(readMM('gse125449/GSE125449_Set1_matrix.mtx'))
barcodes = read.csv("gse125449/GSE125449_Set1_barcodes.tsv",header = FALSE)
features = read.csv("gse125449/GSE125449_Set1_genes.tsv",sep='\t',header = FALSE)

colnames(sc_count) = barcodes[,1]
rownames(sc_count) = features$V2

rownames(meta_ori)=meta_ori[,'Cell.Barcode']
cell_select = meta_ori[meta_ori$Type != 'unclassified','Cell.Barcode']

meta_filt = meta_ori[cell_select,]
table(meta_filt$Type)
sc_count_filt = sc_count[,cell_select]
sc_data = Seurat::CreateSeuratObject(counts = sc_count_filt,
                                     meta.data = meta_filt)
head(sc_data)


######### run RCTD
output_dir = 'result'
output_name = 'gse125449'
cell_type_columns = 'annot'
RCTD_run(sc_data_filt,cell_type_columns,st_count,st_coord,output_dir,output_name)   
SPOTlight_run(sc_data_filt,cell_type_columns,st_count,st_coord,output_dir,output_name)   

for (i in c(1,2,3,4)){
  hcc_st_s = hcc_st_list[[i]]
  st_count = hcc_st_s@assays$RNA@counts
  st_coord = hcc_st_s@meta.data[,c('imagerow','imagecol')]
  output_name = paste('gse125449_st',i,sep="")
  RCTD_run(sc_data,cell_type_columns,st_count,st_coord,output_dir,output_name)   
}

