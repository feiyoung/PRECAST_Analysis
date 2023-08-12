
# Define some functions ---------------------------------------------------

#### Clustering metric
cluster_metric <- function(hy, y, type='ARI'){
  
  require(mclust)
  require(aricode)
  switch(type, 
         ARI= adjustedRandIndex(hy, y),
         NMI = NMI(as.vector(hy), y))
}

# Load data and preprocessing ---------------------------------------------
library(PRECAST)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)
library(dplyr)
name_ID12 <- as.character(c(151507, 151508, 151509, 151510, 151669, 151670,
                            151671, 151672, 151673, 151674, 151675, 151676))


## Read data in an online manner
name_ID <- name_ID12
n_ID <- length(name_ID)
url_sparkA <- "https://github.com/feiyoung/DR-SC.Analysis/raw/main/data/DLPFC_data/brain_"; url_sparkB <- "_spark.Rdata"
url_brainA <- "https://github.com/feiyoung/DR-SC.Analysis/raw/main/data/DLPFC_data/"; url_brainB <- ".rds"
posList <- list()
# genesList <- list()
yList <- list()
spaFeatureList <- list()
seuList <- list()
for(iter in 1:n_ID){
  # iter <- 1
  
  cat('input brain data', iter, '\n')
  # load and read data
  dlpfc <- readRDS(url(paste0(url_brainA, name_ID[iter],url_brainB) ))
  load(url(paste0(url_sparkA, name_ID[iter],url_sparkB) ))
  set.seed(101)
  adjPval <- PvalDF[,2]
  names(adjPval) <- row.names(PvalDF)
  sort_adjPval <- sort(adjPval)
  n <- ncol(dlpfc)
  
  seu1 <- CreateSeuratObject(counts = dlpfc@assays@data$counts, 
                             meta.data = as.data.frame(colData(dlpfc)))
  
  seu1[['RNA']] <- AddMetaData(seu1[['RNA']], metadata=as.data.frame( rowData(dlpfc)))
  seu1@tools$platform <- "Visium"
  seuList[[iter]] <- seu1
  rm(seu1)
  spaFeatureList[[iter]] <- names(sort_adjPval)
  yList[[iter]] <- dlpfc$layer_guess_reordered
  posList[[iter]] <- cbind(dlpfc$row, dlpfc$col) 
}


### Save seuList for 12 data##################
saveRDS(seuList, file='Brain12_seuList.RDS')


### Save gene list for use#######################
selectIntFeatures <- PRECAST:::selectIntFeatures
genelist <- selectIntFeatures(seuList, spaFeatureList)
length(unique(genelist))
write.csv(genelist, file='brain12_genelist2000.csv')

### Save Log-nomalized data list###############
getXList <- function(seuList, genelist){
  require(Seurat)
  set.seed(101)
  seuList <- lapply(seuList, NormalizeData)
  XList <- list()
  nsample <- length(seuList)
  indexList <- list()
  nr <- 0
  yList <- list()
  posList <- list()
  for(i in 1:nsample){
    XList[[i]] <- Matrix::t(seuList[[i]]@assays$RNA@data[genelist,])
    indexList[[i]] <- (nr+1):(nrow(XList[[i]] )+nr)
    nr <- nr + nrow(XList[[i]] )
    y <- as.character(seuList[[i]]$layer_guess_reordered)
    y[is.na(y)] <- "NA"
    yList[[i]] <- y
    posList[[i]] <- cbind(seuList[[i]]$row, seuList[[i]]$col)
  }
  
  return(list(XList=XList, posList=posList, 
              yList=yList, 
              indxList=indexList))
}


datList <- getXList(seuList, genelist)
saveRDS(datList, file='datList.RDS')



# Integration analysis using PRECAST --------------------------------------


sapply(datList$XList, dim)
XList <- datList$XList
posList <- datList$posList
yList <- datList$yList 
indexList <- datList$indxList

K_set <- 2:11; q <- 15
tic <- proc.time() # 
set.seed(1) ## parallel computing
resList <- ICM.EM(XList, posList=posList, q=q, K=K_set, 
                  platform = 'Visium',maxIter = 30,
                  Sigma_equal =F, verbose=T,  coreNum=5, coreNum_int = 5)
toc <- proc.time()
time_used <- toc[3] - tic[3] 


reslist <- SelectModel(resList)

### Calculate the metrics for clustering




ari_idrsc_vec <- sapply(1:12, function(r) cluster_metric(yList[[r]], reslist$cluster[[r]]))
nmi_idrsc_vec <- sapply(1:12, function(r) cluster_metric(yList[[r]], as.vector(reslist$cluster[[r]]), type='NMI'))

anMat <- cbind(ARI=ari_idrsc_vec, NMI=nmi_idrsc_vec)
print(anMat)


ari_idrsc <- cluster_metric(unlist(reslist$cluster), unlist(yList), type='ARI')
nmi_idrsc <- cluster_metric(unlist(reslist$cluster), unlist(yList), type='NMI')
cluster_idrsc <- reslist$cluster
save(ari_idrsc, nmi_idrsc, cluster_idrsc,
     nmi_idrsc_vec, nmi_idrsc_vec, file='Brain12_metric_cluster10_idrsc.rds')




## Get hZ, tSNE and ConCor
hZ_idrsc10 <- matlist2mat(reslist$hZ)
concor_idrsc_cluster10 <- dr_concor(X0, hZ_idrsc10, unlist(yList))

sampleID <- get_sampleID(XList)
set.seed(1)
hZ_tsne_idrsc_cluster10 <- scater::calculateTSNE(t(hZ_idrsc10))
row.names(hZ_tsne_idrsc_cluster10) <- sampleID
save(concor_idrsc_cluster10, hZ_idrsc10, hZ_tsne_idrsc_cluster10,
     file='Brain12_hZtsne_idrsc_cluster10.rds')

save(concor_idrsc_cluster10, file='Brain12_ConCor_idrsc_cluster10.rds')


## calculate hV tSNE for each sample
hVList <- reslist$hV
tsne3_hvlist <- list()
for(r in 1:length(hVList)){
  message("r = ", r)
  tsne3_hvlist[[r]] <- scater::calculateTSNE(t(hVList[[r]]), ncomponents = 3)
}

umap3_hvlist <- list()
for(r in 1:length(hVList)){
  message("r = ", r)
  umap3_hvlist[[r]] <- scater::calculateUMAP(t(hVList[[r]]), ncomponents = 3)
}
save(tsne3_hvlist, umap3_hvlist, file="tsne3_umap3_hvlist_brain12.rds")





# Downstream analysis -----------------------------------------------------


# Combined DEG analysis ---------------------------------------------------

reslist <- selectModel(resList, return_est=T)
housekeep <- read.csv(file='housekeep_genes_brain12.csv')[,2]
Rf <- attr(reslist, 'fit')$Rf
hX <- get_correct_exp(XList, Rf, housekeep)

library(org.Hs.eg.db)
Symbol <- mapIds(org.Hs.eg.db, keys = colnames(hX),
                 keytype = "ENSEMBL", column="SYMBOL")
Symbol[is.na(Symbol)] <- paste0("NA",1: sum(is.na(Symbol)))
colnames(hX) <- Symbol

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

dat_degs <- FindAllMarkers(seuAll)
save(dat_degs, file='housekeep_dat_degs_cluster10.rds')


#Plot DEGs expresion heatMaps based on all intergarted samples
library(dplyr)
n <- 10
dat_degs %>%
  group_by(cluster) %>%
  top_n(n = n, wt = avg_log2FC) -> top10
seu <- seuAll
seu@assays$RNA@var.features <- row.names(seu)
seu[['RNA']]@data[1:4,1:4]
seus <- subset(seu, downsample = 4000)
color_id <- as.numeric(levels(Idents(seus)))
seus <- Seurat::ScaleData(seus)
seus[['RNA']]@scale.data[1:4,1:4]
sum(top10$gene %in%  row.names(seu))

cols_cluster <- c( "#E15759CC", "#59A14FCC", "#FF9DA7CC", "#9C755FCC",
                   "#B07AA1CC", "#4E79A7CC", "#EDC948CC", "#BAB0ACCC",
                   "#ff7f0e", "#aec7e8CC")
features <- top10$gene
## HeatMap
p1 <- doHeatmap(seus, features = features, cell_label= "Domain",
                grp_label = F, grp_color = cols_cluster[color_id],
                pt_size=6,slot = 'scale.data') + 
  theme(legend.text = element_text(size=16),
        legend.title = element_text(size=18, face='bold'),
        axis.text.y = element_text(size=7, face= "italic", family='serif'))
ggsave(paste0('BrainAll',"_top",n,"DEGs_heatmap.pdf"), plot = p1, width = 10, height = 8, units = "in", dpi = 1000)


cutoff <- 0.001 #0.05
dat_degs_sub <- subset(dat_degs,  avg_log2FC >0.25 &   p_val_adj<cutoff )
table(dat_degs_sub$cluster)
sum(table(dat_degs_sub$cluster))

filename <- 'DEGsList_housekeep_jointAnalysis_cutoff0.001.xlsx'
for(k in 1:10){
  # k <- 1
  cat('k=', k, '\n')
  dat_degs_sub3 <- subset(dat_degs,  cluster==k &  p_val_adj<cutoff & avg_log2FC >0.25)
  
  dat <- as.data.frame(dat_degs_sub3)
  dat$gene <- firstup(dat$gene)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Domain', k),
                   append = T, row.names = F)
}



# Combined trajectory inference using Slingshot ------------------------------------

library(slingshot)
load('Brain12_metric_cluster10_idrsc.rds')
load('Brain12_hZtsne_idrsc_cluster10.rds')
cluster_idrsc_vec <- unlist(cluster_idrsc)
sds <- slingshot(hZ_idrsc10, cluster_idrsc_vec,start.clus = 5) 

save(sds, file = 'Brain12_traject_sds.rds')

pt <- slingPseudotime(sds)
common.pseudo <- rowMeans(pt, na.rm=TRUE)


## Test the pseudotime-associated genes###############
### Pearson correlation test###################
save(hX, file='batch_corrected_exprs_iDRSC.rds')
corrTestList <- list()
for(j in 1:2000){
  # j <- 1
  message("j = ",j)
  corrTestList[[j]] <-  cor.test(hX[,j], common.pseudo)
}
names(corrTestList) <- colnames(hX)

save(corrTestList, file='cor_expre_pseudotime_brain12.rds')
corrTestMat <- matrix(NA, 2000, 4)
colnames(corrTestMat) <- c("pearson_corr", "t statistic", "p_val", "adjp_val")
row.names(corrTestMat) <- colnames(hX)
for(j in 1:2000){
  message("j = ",j)
  corTest <- corrTestList[[j]]
  corrTestMat[j,1:3] <- c(corTest$estimate, corTest$statistic, corTest$p.value)
}
corrTestMat[, 4] <- p.adjust(corrTestMat[,3],method="bonferroni")
head(corrTestMat)
sum(corrTestMat[,4]<0.05 & abs(corrTestMat[,1])>0.1)
sort(corrTestMat[,1], decreasing = T)[1:20]
idx <- which(corrTestMat[,4]<0.05 & abs(corrTestMat[,1])>0.1)
write.csv(corrTestMat[idx, ], file='corrTestMat_sig.csv')

corrSig <- read.csv(file='corrTestMat_sig.csv', row.names = 1)
head(corrSig)
cutoff <- 0.001

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
row.names(corrSig) <- firstup(row.names(corrSig))

corrSig_new <- corrSig[corrSig[,4]<cutoff,]
write.csv(corrSig_new, file='corrTestMat_sig_001_brain12.csv')

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
summary(df1$P.DE)
summary(df1$nlog10P)
df1$ID <- 1:nrow(df1)
library(ggplot2)
p = ggplot(df1,aes(ID, nlog10P, color=Ont))+
  geom_point(aes(size=DE, color=Ont), alpha=0.5) +labs(color='Category', 
                                                       x="Pathway ID",y=expression(-log[10](P)),title=NULL)

cols <- c("steelblue3", "goldenrod", "brown3")
label <- rbind(subset(df1, Ont=='BP')[id_label_bp[c(2,3)],],
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



# Conditional SVA analysis ------------------------------------------------

seulist <- readRDS(file='Brain12_seulist.RDS')
indexList <- get_indexList(posList)

conditionalSVGs <- function(seu, hZ, genelist=NULL, num_core=10, verbose=T){
  require(SPARK)
  require(Seurat)
  require(Matrix)
  require(doParallel)
  
  if(!is.null(genelist)) seu <- seu[genelist, ]
  set.seed(101)
  mat <- seu[["RNA"]]@counts
  pos <- cbind.data.frame(seu$row, seu$col)
  rownames(pos) <- colnames(mat)
  ## filter genes and cells/spots and
  spark_brain <- CreateSPARKObject(counts=mat,
                                   location= pos,
                                   percentage = 0,
                                   min_total_counts = 0)
  
  
  spark_brain@lib_size <- apply(spark_brain@counts, 2, sum)
  tic <- proc.time()
  spark_brain <- spark.vc(spark_brain, covariates = hZ, lib_size = spark_brain@lib_size,
                          num_core = num_core,  fit.model = 'gaussian', verbose = verbose)
  toc1 <- proc.time() - tic
  
  tic <- proc.time()
  spark_brain <- spark.test(spark_brain, check_positive = T, verbose = verbose)
  toc2 <- proc.time() - tic
  
  PvalDF <- spark_brain@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
  
  attr(PvalDF, 'time_used') <- toc2 + toc1
  return(PvalDF)
} 

## SVGs test for each data ##########################
SVGsList <- list()
for(r in 1:12){
  res_tmp <- conditionalSVGs(seulist[[r]], hZ_idrsc10[indexList[[r]], ])
  SVGsList[[r]] <- res_tmp
  save(SVGsList, file='brain12_SVGs_SVGsList.rds') # use this method to savd data in time
}

### write the SVGs into xlsx file ########################
library(org.Hs.eg.db)
filename <- 'Conditional_SVGs_brain12.xlsx'
for(r in 1:12){
  #  r <- 1
  cat('r=', r, '\n')
  svgs_tmp <- subset(SVGsList[[r]],  adjusted_pvalue<0.01)
  Symbol <- mapIds(org.Hs.eg.db, keys = row.names(svgs_tmp), 
                   keytype = "ENSEMBL", column="SYMBOL")
  svgs_tmp$symbol <- Symbol
  dat <- as.data.frame(svgs_tmp)
  dat$symbol <- firstup(dat$symbol)
  xlsx::write.xlsx(dat, file=filename, sheetName = paste0('Brain', r),
                   append = T)
}


### SVGs enrichment analysis
library(gprofiler2)
#r <- 1
termList <- list()
for(r in 1:12){
  cat("r = ", r, '\n')
  que1 <- sigGeneList[[r]]
  gostres <- gost(query = que1,
                  organism = "hsapiens", correction_method='fdr')
  termList[[r]] <- gostres
}
save(termList, file='profiler_termList_brain12.rds')



## Bar plot---------------------------
get_top_pathway1 <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id","p_value")])
  }
  
  return(df_sub)
}
barPlot_enrich <- function(top_dat, source='Ont', term_name="Term", nlog10P='nlog10P',
                           bar_width=0.8, base_size=20, font_family='serif', cols= ggthemes::canva_pal()(4)){
  # source='Ont'; term_name="Term"; nlog10P='nlog10P'
  require(ggplot2) # y=term_name,
  order_idx <- order(top_dat[,nlog10P])
  top_dat <- top_dat[order_idx,]
  top_dat[, term_name] <- factor(top_dat[, term_name], levels=top_dat[order_idx,term_name])
  p1 <- ggplot(data=top_dat, aes_string(x=term_name,y=nlog10P, fill=source)) +
    scale_fill_manual(values=cols)+
    geom_bar(position = "dodge", stat="identity",width =bar_width)+ coord_flip() +
    theme_classic() + theme(text=element_text(size=base_size, family=font_family)) 
  return(p1)
}

dat1$nlog10P <- -log10(dat1$p_value)
barPlot_enrich(subset(dat1, source=='KEGG'), source='source', 'term_name', 'nlog10P')
source_set <- c("GO:BP","GO:CC", "GO:MF",   'KEGG', "HPA")
## KEGG
source1 <- source_set[3]
ss <- which(source_set==source1)
ntop = 10
cols <- c("steelblue3", "goldenrod", "brown3", "#f98866", "#CE6DBD")
names(cols) <- source_set
pList_enrich <- list()
for(ii in 1:12){
  ## ii <- 1
  message("ii=", ii)
  gostres2 <- termList[[ii]]
  dat_tmp <- gostres2$result
  table(dat_tmp$source)
  dat1 <- get_top_pathway1(dat_tmp[dat_tmp$term_id!= "KEGG:05171",], ntop=ntop)
  dat1$nlog10P <- -log10(dat1$p_value)
  dat1 <- subset(dat1, source==source1)
  pList_enrich[[ii]] <- barPlot_enrich(dat1[order(dat1$nlog10P),], source='source', 'term_name',
                                       'nlog10P', cols=cols[source_set[ss]]) + ylab("-log10(p-adj)")+
    xlab("Biological terms") + ggtitle(paste0("Sample", ii)) + theme(legend.position = 'none')
  
}

# p12 <- cowplot::plot_grid(plotlist = pList_enrich, nrow=6, ncol=2)
p12 <- patchwork::wrap_plots(plotlist = pList_enrich, nrow=6, ncol=2)
if(source1%in% source_set[1:3]){
  
  ggsave(file=paste0("Brain12_SVG_enrich", strsplit(source1, split=':')[[1]][2], "_barplot.png"), plot = p12,
         width = 25, height =24, units = "in", dpi = 400)
  
}else{
  ggsave(file=paste0("Brain12_SVG_enrich", source1, "_barplot.png"), plot = p12,
         width = 25, height =24, units = "in", dpi = 400)
}
