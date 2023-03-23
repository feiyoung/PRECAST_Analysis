

# Fig 2a ------------------------------------------------------------------


r <- 10
pos10 <- posList[[r]]

for(i in 1:4){
  message("i = ", i)
  pos10 <- cbind(pos10,hZ_umap3List_allsample[[i]][indexList[[r]], ])
}

dim(pos10)
head(pos10)


for(j in 1:4){
  cluster_tmp <- clusterMat_sub[indexList[[r]],j]
  pos10 <- cbind(pos10,cluster_tmp)
}
colnas <- c("coord x", "coord y", rep(names(hZ_umap3List_allsample[1:4]), each=3), names(hZ_umap3List_allsample[1:4]))
colnames(pos10) <- colnas
getwd()
write.csv(pos10, file='tmp.csv')

##### 2b #####
dat_tmp <- Reduce(cbind,c(tsne2List[1:4], tsne2List[9]))
dat_tmp <- cbind(dat_tmp, meta_data[,c(1,3)])
colnas <- c(rep(names(tsne2List[c(1:4, 9)]), each=2), "Sample", "Domain")
colnames(dat_tmp) <- colnas
write.csv(dat_tmp, file='tmp.csv', row.names = F)

##### 2c #####
ariMat_sub <- cbind(ariMat_sub, "PASTE"=ari_paste)
row.names(ariMat_sub) <- paste0("Sample", 1:12)
write.csv(ariMat_sub, file='tmp.csv')

##### 2d #####
write.csv(corMat, file='tmp.csv')

##### 2e&f #####

dat <- as.data.frame(seu[["Spatial"]]@cell.embeddings)
genes <- t(seu[['RNA']]@scale.data[union(toupper(gene_eachdomain), toupper(gene_assoc_pseudotime)),])
dat_tmp <- cbind(dat, genes)
write.csv(dat_tmp, file='tmp.csv')


##### 2g #####
write.csv(df1, file='tmp.csv')

##### 2h #####
dat_tmp <- gostres2$result
write.csv(as.matrix(dat_tmp), file='tmp.csv')

##### 3a #####

ARI_tmp <- rbind(ariMat[,main_order_names2], ariVec[main_order_names2])
row.names(ARI_tmp) <- c(paste0("Sample", 1:8), "Combined")
write.csv(ARI_tmp,  file='tmp.csv')
cLISI_tmp <- clisiMat[,main_order_names2]

iLISI_tmp <- ilisiMat[,main_order_names2]

write.csv(cbind(cLISI_tmp, iLISI_tmp),  file='tmp.csv')


##### 3b #####
dat_tmp <- Reduce(cbind,c(hZtsneList[c(1,3,4, 8)], hZtsneList[10]))
dat_tmp <- cbind(dat_tmp, meta_data[,c(1,2)])
colnas <- c(rep(names(hZtsneList[c(1,3,4, 8, 10)]), each=2), "Sample", "Domain")
colnames(dat_tmp) <- colnas
write.csv(dat_tmp, file='tmp.csv', row.names = F)

##### 3c #####
dat <- as.data.frame(t(seu[["RNA"]]@scale.data))
dat$Domain <- Idents(seu)
write.csv(dat, file='tmp.csv')


##### 3d #####
Mac_tmp <- round(MacR2Vec[main_order_names2], 2)
write.csv(Mac_tmp, file='tmp.csv')


##### 3e #####
datList_here <- list()
for(slice in 1:3){
  
  for(jj in 1:6){
    if(jj==1){
      dats <- dat_List_all[[slice]][[jj]][1:3]
    }else{
      tmp <- dat_List_all[[slice]][[jj]][3]
      dats <- cbind(dats, tmp)
    }
   
  }
  datList_here[[slice]] <- dats
}

dat_all <- datList_here[[1]]
dat_all$sample <- 1
for(r in 2:3){
  tmpdat <- datList_here[[r]]
  tmpdat$sample <- r
  dat_all <- rbind(dat_all, tmpdat)
}
write.csv(dat_all, file='tmp.csv')


##### 3f #####
dat_tmp <- dat1[,c(1:2,4)]
write.csv(dat_tmp, file='tmp.csv', row.names = F)
##### 3g #####

dat_tmp <- as.data.frame(t(logcounts(sce.liver)[features_reorder,]))

dat_tmp$Pseudotime  <- range01(sce.liver$Pseudotime)
dat_tmp$Domain <- sce.liver$domain 

write.csv(dat_tmp, file='tmp.csv')


##### 4b #####
head(datList[[1]])
dat_tmp <- as.data.frame(datList[[1]])
colnames(dat_tmp) <- c("Coord y", "Coord y", "Domain")
dat_tmp$Sample <- 1
for(r in 2:16){
  message("r = ", r)
  dat_tmp2 <- as.data.frame(datList[[r]])
  colnames(dat_tmp2) <- c("Coord y", "Coord y", "Domain")
  dat_tmp2$Sample <- r
  dat_tmp <- rbind(dat_tmp, dat_tmp2)
}
write.csv(dat_tmp, file='tmp.csv')

##### 4c #####

dat_tmp <- Reduce(cbind,c(tsne2List[1:4], tsne2List[9]))
dat_tmp <- cbind(dat_tmp, meta_data[,c(1,2)])
colnas <- c(rep(names(tsne2List[c(1:4, 9)]), each=2), "Sample", "Domain")
colnames(dat_tmp) <- colnas
write.csv(dat_tmp, file='tmp.csv', row.names = F)

##### 4d #####
dat_tmp <- percentage_long
write.csv(dat_tmp, file='tmp.csv', row.names = F)

##### 4e #####
dat_tmp <-  rbind(MacR2Vec[main_order_names2], ariMat[,main_order_names2])
write.csv(dat_tmp, file='tmp.csv', row.names = F)

##### 4f #####
dat_tmp <- Reduce(rbind, datList[1:8])
write.csv(dat_tmp, file='tmp.csv')



##### 5b #####
dat_tmp <- Reduce(rbind, datList)
write.csv(dat_tmp, file='tmp.csv')

##### 5c #####
dat_tmp <- Reduce(cbind,c(hZtsneList[1:4], hZtsneList[9]))
dat_tmp <- cbind(dat_tmp, meta_data[,c(1,2)])
colnas <- c(rep(names(hZtsneList[c(1:4, 9)]), each=2), "Sample", "Domain")
colnames(dat_tmp) <- colnas
write.csv(dat_tmp, file='tmp.csv', row.names = F)


#### 5d #########
for(r in 1:4){
  datList[[r]][,"Proportion2"] <- datList0[[r]]$Proportion
  datList[[r]][,"Proportion3"] <- datList2[[r]]$Proportion
}
for(r in 1:4){
  datList[[r]][,"Sample"] <- paste0("HCC", r)
  
}
dat_tmp <- Reduce(rbind, datList)
write.csv(dat_tmp, file='tmp.csv')

##### 5e #####
dat_tmp <- Reduce(rbind, datList)
write.csv(dat_tmp, file='tmp.csv')
