iDR_SC_newname <- "PRECAST"


##############color related settings
colorbar_adj_transparent <- function(colors, alpha=0.6, plot=T){
  require(colorspace)
  ramp.list = adjust_transparency(colors,   alpha = alpha)
  print(ramp.list)
  if(plot==T){
    barplot(rep(1, length(ramp.list)), axes = FALSE, space = 0, col = ramp.list)
  }
  return(ramp.list)
}

# cols_alpha <- pal
cols_alpha <- c("#E04D50", "#4374A5", "#F08A21","#2AB673", "#FCDDDE",  "#70B5B0", "#DFE0EE" ,"#D0B14C")


cols_cluster <- c(cols_alpha[1:2], rep(cols_alpha[-c(1:2)], each=2))
cols_cluster[seq(4, 14, by=2)] <- colorbar_adj_transparent(cols_cluster[seq(4, 14, by=2)], alpha=0.5)

cols_alpha2 <-  c("#E04D50", "#4374A5", "#F08A21","#2AB673", "#FCDDDE",  "#70B5B0", "#DFE0EE" ,"#D0B14C", '#57A6D9')
cols_cluster2 <- c(cols_alpha2[1:2], rep(cols_alpha2[-c(1:2)], each=2))
cols_cluster2[seq(4, 16, by=2)] <- colorbar_adj_transparent(cols_cluster2[seq(4, 16, by=2)], alpha=0.5)



## set the order of methods
main_order_names <- c("iDR-SC", "Seurat V3",  "Harmony","fastMNN",
                      "Scanorama", "scGen","scVI","MEFISTO", "PASTE")
#cluster order 
order_names <- c("iDR-SC", "Seurat V3" ,"Harm-SC-MEB", "Harm-Louvain", 
                 "fastMNN-SC-MEB", "fastMNN-Louvain", "Scanorama-SC-MEB", "Scanorama-Louvain",
                 "scGen-SC-MEB", "scGen-Louvain","scVI-SC-MEB", "scVI-Louvain",
				 'MEFISTO-SC-MEB', 'MEFISTO-Louvain', 'PASTE-SC-MEB', 'PASTE-Louvain')
					  
## select colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
get_trans_colors <- function(color, num=2){
  require(colorspace)
  alphavec <- seq(0.2, 1, length=num)
  color_vec <- rep(NA, num)
  for(i in 1:num)
    color_vec[i] = adjust_transparency(color,   alpha = alphavec[i])
  return(color_vec)
}

# 
color2bar_gradient <- function(color1, color2, len_grad=10, plot=T){
  require(grDevices)
  ramp <- colorRamp(c(color1, color2))
  ramp.list <- rgb( ramp(seq(0, 1, length = len_grad)), max = 255)
  if(plot==T){
    barplot(rep(1, length(ramp.list)), axes = FALSE, space = 0, col = ramp.list)
  }
  print(ramp.list)
  return(ramp.list)
  
}
tune_colors <- function(colors, alpha=0.6, plot=T){
  require(colorspace)
  ramp.list = adjust_transparency(colors,   alpha = alpha)
  
  if(plot==T){
    barplot(rep(1, length(ramp.list)), axes = FALSE, space = 0, col = ramp.list)
  }
  return(ramp.list)
}



################# tool functions

range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

firstup <- function(x) {
## First letter use upper capital
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

get_sampleID <- function(XList){
  sampleID <- list()
  r_max <- length(XList)
  for(r in 1:r_max){
    sampleID[[r]] <- rep(r, nrow(XList[[r]]))
  }
  sampleID <- unlist(sampleID)
  return(sampleID)
}


get_indexList <- function(alist){
  nsample <- length(alist)
  nr <- 0
  indexList <- list()
  for(i in 1:nsample){
    indexList[[i]] <- (nr+1):(nrow(alist[[i]] )+nr)
    nr <- nr + nrow(alist[[i]] )
  }
  return(indexList)
}

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


replace_cluster <- function(clusters, new_names){
  
  uni_clus <- unique(clusters)
  new_names <- new_names[uni_clus]
  n_clus <- length(unique(clusters)) 
  new_vec <- clusters
  for(i in 1:n_clus){
    new_vec[clusters==uni_clus[i]] <- new_names[i]
  }
  return(new_vec)
}



get_top_pathway1 <- function(df, ntop=10, source_set=c("GO:MF", "GO:CC", "GO:BP", 'KEGG', "HPA")){
  df_sub <- NULL
  n_source <- length(source_set)
  for(i in 1:n_source){
    tmp <- subset(df, source==source_set[i])
    df_sub <- rbind(df_sub, tmp[1:ntop, c("source", "term_name", "term_id","p_value")])
  }
  
  return(df_sub)
}

coordinate_rotate <- function(pos, theta=0){# counter-clock rotation
  pos_new <- pos
  pos_new[,1] <- pos[,1]*cos(theta) - pos[,2]*sin(theta)
  pos_new[,2] <- pos[,1]*sin(theta) + pos[,2]*cos(theta)
  return(pos_new)
}  
rotate_angles <- c(90.1, -40.5)


align_colors <- function(ClusterMat,base_cluster, cols_cluster){
  # j <- 1
  if(length(cols_cluster) < max(apply(ClusterMat, 2, max))) stop("Number of cols is not enough!")
  givecolors <- function(mapID, base_cols_cluster, other_cols){
    m <- ncol(mapID)
    color_ID <- unname(base_cols_cluster)
    row_clusterID <- as.numeric(names(mapID))
    row_clusterCols <- rep("", m)
    col_used_count <- rep(0, length(base_cols_cluster))
    tmp_base_cluster <- NULL
    for(i in 1: m){
      # i <- 1
      
      if((mapID[1,i] %in% color_ID) &&  (max(mapID[2, mapID[1,] == mapID[1,i]]) == mapID[2,i])){
        col_chose <- names(base_cols_cluster)[mapID[1, i]]
        row_clusterCols[i] <- col_chose
        tmp_base_cluster <- c(tmp_base_cluster, mapID[1, i])
        color_ID <- setdiff(color_ID, mapID[1, i])
        
      }else{
        row_clusterCols[i] <- other_cols[1]
        other_cols <- other_cols[-1]
      }
      
      col_used_count[mapID[1, i]] <- col_used_count[mapID[1, i]] + 1
      
      
    }
    if(any(is.na(row_clusterCols))) row_clusterCols[is.na(row_clusterCols)] <- names(base_cols_cluster)[color_ID]
    return(row_clusterCols)
  }
  base1 <- sort(unique(base_cluster))
  n_base <- length(base1)
  base_cols_cluster <- base1
  names(base_cols_cluster) <- cols_cluster[1:n_base]
  other_cols <- cols_cluster[-(1:n_base)]
  colorList <- list()
  for(j in 1:ncol(ClusterMat)){
    # j <- 2
    stab <- table(ClusterMat[,j], base_cluster)
    mapID <- rbind(apply(stab, 1, which.max), apply(stab, 1, max))
    colorList[[j]] <- givecolors(mapID, base_cols_cluster, other_cols)
  }
  return(colorList)
}

subsample <- function(hZ_idrsc, cluster_idrsc_vec, sample_rate = 0.5){
  nn <- length(cluster_idrsc_vec)
  idx <- sort(sample(nn, nn* sample_rate))
  return(list(hZ=hZ_idrsc[idx, ], cluster=cluster_idrsc_vec[idx]))
}

# Plot functions ----------------------------------------------------------
pt_size_sample <- 0.2; pt_alpha_sample <- 0.7
pt_size_cluster <- 0.6; pt_alpha_cluster <- 0.5
base_axis_size <- 20
width_tsne_fig_nolengend <- 6.6; height_tsne_fig_nolengend <- 6
bar_plot_hip2 <- function(percentage_long,geom_bar_position,legend_position, fill_color,
                     title_name="Hippo"){
  ggplot(percentage_long, aes(y = value, x = factor(Cluster), fill = variable)) +        ## global aes
    scale_fill_manual(values= fill_color,name = 'Cell Type')+
    # ggthemes::tableau_color_pal("Tableau 20")(ncol(percentage)-1)
    geom_bar(position=geom_bar_position, stat="identity",width=0.7,color="black") +
    ggtitle(title_name) +
    theme_bw()+xlab("")+ylab("")+
    theme_classic() +
    theme(plot.title = element_text(size = 14,hjust = 0.5),
          text = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "grey"),
          legend.position = legend_position,
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "gray13"))
}


#' @title Visualize RGB plot from tSNE.
#' @description We summarized the inferred low dimensional components into three tSNE components and visualized the three resulting components with red/green/blue (RGB) colors in the RGB plot.
#' @param location A n by k location matrix. n is spot number.
#' @param latent_dat A d by n matrix of low dimensional components.
#' @param pointsize The point size of each spot.
#' @param textsize The text size in the legend.
#' @return A list.
#' \item{RGB}{A data frame with five columns: x coordinate, y coordinate, R, G, and B color index}
#' \item{figure}{A ggplot object for RGB plot from tSNE}
#'
#' @import Rtsne
#'
#' @export
plot_RGB=function(location, embed_3d, pointsize=2,textsize=15){

  suppressMessages(require(ggplot2))

  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")


  r = (embed_3d[,1]-min(embed_3d[,1]))/(max(embed_3d[,1])-min(embed_3d[,1]))
  g = (embed_3d[,2]-min(embed_3d[,2]))/(max(embed_3d[,2])-min(embed_3d[,2]))
  b = (embed_3d[,3]-min(embed_3d[,3]))/(max(embed_3d[,3])-min(embed_3d[,3]))
  x =  info$sdimx
  y =  info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
      geom_point(size=pointsize) +
      scale_color_identity()+
      theme_void()+
      theme(plot.title = element_text(size = textsize),
              text = element_text(size = textsize),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 22) ,
              legend.position = "bottom")

  p1
}


rotate_90_clockwise <- function(pl){
  require(ggplot2)
  pl + coord_flip() + scale_x_reverse()
}

plot_RGB_tSNE <- function (location, latent_dat, pointsize = 2, textsize = 15) 
{
    info = as.data.frame(location)
    colnames(info) = c("sdimx", "sdimy")
    PCvalues = latent_dat
    tsne <- Rtsne(t(PCvalues), dims = 3, check_duplicates = FALSE)
    r = (tsne$Y[, 1] - min(tsne$Y[, 1]))/(max(tsne$Y[, 1]) - 
        min(tsne$Y[, 1]))
    g = (tsne$Y[, 2] - min(tsne$Y[, 2]))/(max(tsne$Y[, 2]) - 
        min(tsne$Y[, 2]))
    b = (tsne$Y[, 3] - min(tsne$Y[, 3]))/(max(tsne$Y[, 3]) - 
        min(tsne$Y[, 3]))
    x = info$sdimx
    y = info$sdimy
    dat = data.frame(x, y, r, g, b)
    p1 = ggplot(data = dat, aes(x = x, y = y, col = rgb(r, g, 
        b))) + geom_point(size = pointsize) + scale_color_identity() + 
        ggtitle(paste0("RGB tSNE")) + theme_void() + theme(plot.title = element_text(size = textsize), 
        text = element_text(size = textsize), legend.position = "bottom")
    return(list(RGB = dat, figure = p1))
}

plot_RGB_UMAP <- function (location, latent_dat, pointsize = 2, textsize = 15) 
{
    require(umap)
    info = as.data.frame(location)
    colnames(info) = c("sdimx", "sdimy")
    PCvalues = latent_dat
    umap <- umap(t(PCvalues), n_components = 3)
    r = (umap$layout[, 1] - min(umap$layout[, 1]))/(max(umap$layout[, 
        1]) - min(umap$layout[, 1]))
    g = (umap$layout[, 2] - min(umap$layout[, 2]))/(max(umap$layout[, 
        2]) - min(umap$layout[, 2]))
    b = (umap$layout[, 3] - min(umap$layout[, 3]))/(max(umap$layout[, 
        3]) - min(umap$layout[, 3]))
    x = info$sdimx
    y = info$sdimy
    dat = data.frame(x, y, r, g, b)
    p1 = ggplot(data = dat, aes(x = x, y = y, col = rgb(r, g, 
        b))) + geom_point(size = pointsize) + scale_color_identity() + 
        ggtitle(paste0("RGB UMAP")) + theme_void() + theme(plot.title = element_text(size = textsize), 
        text = element_text(size = textsize), legend.position = "bottom")
    return(list(RGB = dat, figure = p1))
}


# barPlot_real(lisiVec[main_order_names], ylabel = 'cLISI', cols=cols_alpha)+coord_polar()
# 环形柱状图

barPlot_real <- function(vec, ylabel='ARI', cols=NULL,...){
  require(ggplot2)
  
  
  ## filter vec
  N <- length(vec)
  vec_use <-vec[!is.na(vec)]
  
  
  df_use <- data.frame(value=vec_use, 
                       Method=names(vec_use))
  df_use$Method <- factor(df_use$Method, levels=names(vec_use))
  
  
  
  
  ## CCor
  p1 <- ggplot(df_use, aes(x=Method, y=value, fill=Method)) + 
    geom_bar(position = "dodge", stat="identity",width = 1, ...) + # , ...
    #geom_errorbar( aes(ymin=value-sd, ymax=value+sd), width=0.4, colour="orange",  size=1.3, position=position_dodge(.9)) + 
    #facet_grid(beta~Error , scales="fixed",labeller = label_bquote(beta == .(beta))) 
    labs(y=ylabel, x=NULL)+ 
    scale_x_discrete(breaks = NULL) 
  
  if(is.null(cols)){
    return(p1)
  }else{
    p1 + scale_fill_manual(values = cols)
  }
}


plot_scatter <- function (
    embed_use, meta_data, label_name, xy_names=c('tSNE1', 'tSNE2'), no_guides = FALSE, 
    palette_use = tableau_color_pal()(10), 
    pt_size = 4, point_size = 0.5, point_alpha=1, 
    base_size = 12, do_points = TRUE, do_density = FALSE, border_col='gray',
    legend_pos='right', legend_dir='vertical') {
  require(dplyr)
  require(ggthemes)
  require(ggrepel)
  require(data.table)
  plt_df <- embed_use %>% data.frame() %>% cbind(meta_data) %>% 
    dplyr::sample_frac(1L)
  plt_df$given_name <- plt_df[[label_name]]
  
  plt <- plt_df %>% ggplot(aes_string(colnames(plt_df)[1],colnames(plt_df)[2], col = label_name, 
                                      fill = label_name)) + #  + theme_tufte(base_size = base_size, ticks= show_ticks)
    theme(axis.text.x=element_text(size=base_size, color=1),
          axis.text.y=element_text(size=base_size, color=1),
          axis.title.x = element_text(size=base_size+2, color='black'),
          axis.title.y = element_text(size=base_size+2, color='black'),
          strip.text =  element_text(size=base_size, color='black'),
          strip.background = element_rect(
            linetype = 'solid', color='gray3'
          ),
          legend.direction = legend_dir, legend.position = legend_pos,
          legend.text=element_text(size=base_size+1),
          legend.title=element_text(size=base_size+2),
          panel.background= element_rect(fill = 'white', color=border_col))+
    guides(color = guide_legend(override.aes = list(stroke = 1, 
                                                    alpha = 1, shape = 16, size = 4)), alpha = FALSE) + 
    scale_color_manual(values = palette_use) + scale_fill_manual(values = palette_use) + 
    theme(plot.title = element_text(hjust = 0.5)) + labs(x = xy_names[1], 
                                                         y = xy_names[2])
  if (do_points) 
    plt <- plt + geom_point( size = point_size, alpha=point_alpha)
  if (do_density) 
    plt <- plt + geom_density_2d()
  if (no_guides) 
    plt <- plt + guides(col = 'none', fill = 'none', alpha = 'none')
  
  return(plt)
}



doHeatmap <- function(seu, features=NULL, cell_label='Cell type', grp_label = FALSE,
                      pt_size=4, grp_color=NULL, ...){
  require(ggplot2)
  ngrp <- nlevels(Idents(seu))
  if(is.null(grp_color)){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    grp_color <- gg_color_hue(ngrp)
  }
  
  
  Seurat::DoHeatmap(object = seu, features=features, group.colors = grp_color[1:ngrp], label = grp_label, ...) +
    guides(color = guide_legend(title=cell_label,override.aes = list(stroke = 1, 
                  alpha = 1, shape = 16, size = pt_size, color=grp_color[1:ngrp])), 
           alpha =  "none")
}



boxPlot_real <- function(mat, ylabel='ARI', cols=NULL, ...){
  
  require(ggplot2)
  ## filter mat
  N <- nrow(mat)
  mat_use <- mat[,which(colSums(is.na(mat)) != N)]
  
  
  df_use <- data.frame(value=as.vector(mat_use), Method=rep(colnames(mat_use), each=N))
  df_use$Method <- factor(df_use$Method, levels=colnames(mat_use))
  p1 <- ggplot(df_use, aes(x=Method, y=value, fill=Method)) +
    geom_boxplot(...) +  
    # facet_grid(beta~Error,scales= "fixed",labeller = label_bquote(beta == .(beta)) )+ 
    labs(x=NULL, y= ylabel) + scale_x_discrete(breaks = NULL)
  if(is.null(cols)){
    return(p1)
  }else{
    p1 + scale_fill_manual(values = cols)
  }
  
  
}

volinPlot_real <- function(mat, ylabel='ARI', cols=NULL,...){
  require(ggplot2)
  ## filter mat
  N <- nrow(mat)
  mat_use <- mat[,which(colSums(is.na(mat)) != N)]
  
  
  df_use <- data.frame(value=as.vector(mat_use), Method=rep(colnames(mat_use), each=N))
  df_use$Method <- factor(df_use$Method, levels=colnames(mat_use))
  
  p1 <- ggplot(df_use, aes(x = Method, y = value, fill = Method)) + 
    geom_violin(aes(fill = Method ), color = "transparent", alpha = 0.5) +
    geom_boxplot(outlier.alpha = 0, coef = 0, color = "gray40", width = 0.4) +
    labs(x = "", y = ylabel) +
    theme_classic() + theme( axis.text.x = element_blank() ) + theme(text = element_text(size=20))
  if(is.null(cols)){
    p1
  }else{
    pal <- cols
    p1+scale_fill_manual(values = pal, name = "") 
  }
 
}


barPlot_real2 <- function(vec, ylabel='ARI',cols=NULL,line_size=1, pt_size=2){
  require(ggplot2)
  
  ## filter vec
  N <- length(vec)
  vec_use <-vec[!is.na(vec)]
  
  
  df_use <- data.frame(value=vec_use, 
                       Method=names(vec_use))
  df_use$Method <- factor(df_use$Method, levels=names(vec_use))
  
  p1 <- ggplot(df_use,
               aes(x = value,
                   xend = 0,
                   y = Method,
                   yend = Method,
                   colour = Method)) +
    geom_segment(size=line_size) +
    geom_point(size=pt_size) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_discrete(limits = rev) + coord_flip()+
    labs(x=ylabel, y=NULL)+theme_classic() + theme( axis.text.x = element_blank() ) + 
    theme(text = element_text(size=20))
  if(is.null(cols)){
    return(p1)
  }else{
    p1 + scale_fill_manual(values = cols)
  }
}




do_scatter <- function (
    umap_use, meta_data, label_name, no_guides = FALSE, 
    do_labels = FALSE, nice_names, palette_use = tableau_color_pal()(10), 
    pt_size = 4, point_size = 0.5, pt_shape = '.',
    base_size = 12, do_points = TRUE, do_density = FALSE) {
    plt_df <- umap_use %>% data.frame() %>% cbind(meta_data) %>% 
        dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]
    if (!missing(nice_names)) {
        plt_df %<>% dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))
        plt_df[[label_name]] <- plt_df$nice_name
    }
    plt <- plt_df %>% ggplot(aes_string("X1", "X2", col = label_name, 
        fill = label_name)) + theme_tufte(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, 
            alpha = 1, shape = 16, size = 4)), alpha = 'none') + 
        scale_color_manual(values = palette_use) + scale_fill_manual(values = palette_use) + 
        theme(plot.title = element_text(hjust = 0.5)) + labs(x = "UMAP 1", 
        y = "UMAP 2")
    if (do_points) 
        plt <- plt + geom_point(shape = pt_shape, size = point_size)
    if (do_density) 
        plt <- plt + geom_density_2d()
    if (no_guides) 
        plt <- plt + guides(col = 'none', fill = 'none', alpha = 'none')
    if (do_labels) {
        plt <- plt + geom_text_repel(data = data.table(plt_df)[, 
            .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
            label.size = NA, aes_string(label = label_name), 
            color = "black", size = pt_size, alpha = 1, segment.size = 0) + 
            guides(col = 'none', fill = 'none')
    }
    return(plt)
}

mytheme_graybox <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                             base_rect_size = base_size/22, border_color = 'gray10')
{
  half_line <- base_size/2
  t <- theme(panel.background = element_rect(fill = "white",
                                             colour = NA), panel.border = element_rect(fill = NA,
                                                                                       colour = border_color),
             #line = element_blank(), #rect = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE,
             panel.spacing = unit(half_line/2, "pt"), plot.margin = unit(rep(0.2,4), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t
}


mytheme_void <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                          base_rect_size = base_size/22)
{
  half_line <- base_size/2
  t <- theme(panel.background = element_rect(fill = "white",
                                             colour = NA), panel.border = element_rect(fill = NA,
                                                                                       colour = "grey"),
             line = element_blank(), #rect = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE, panel.spacing = unit(half_line,
                                                                                                                                                                           "pt"), plot.margin = unit(c(0, 0, 0, 0), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t
}


mytheme <- function(legend.direction = "horizontal",
                    legend.position = "bottom", type = 'rect'){
  
  if(type=='rect'){
    th <- theme(axis.text.x=element_text(size=16, color=1, face='plain'),
                axis.text.y=element_text(size=16, color=1, face='plain'),
                axis.title.x = element_text(size=18, color='black', face='plain'),
                axis.title.y = element_text(size=18, color='black',face='plain'),
                strip.text =  element_text(size=16, color='black', face='plain'),
                strip.background = element_rect(
                  linetype = 'solid', color='gray3'
                ),
                legend.direction = legend.direction, legend.position = legend.position,
                legend.text=element_text(size=17, face='plain'),
                legend.title=element_text(size=18, face='bold'),
                panel.background= element_rect(fill = 'white', color='gray'))
  }
  return(th)
}



doHeatmap.matrix <- function(corMat, cluster_orderd,legend_title='Cell type', grp_label = FALSE,
                             pt_size=4, grp_color=NULL){
  
  doHeatmap <- function(seu, features=NULL, cell_label='Cell type', grp_label = FALSE,
                        pt_size=4, grp_color=NULL, ...){
    require(ggplot2)
    ngrp <- nlevels(Idents(seu))
    if(is.null(grp_color)){
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      grp_color <- gg_color_hue(ngrp)
    }
    
    
    Seurat::DoHeatmap(object = seu, features=features, group.colors = grp_color[1:ngrp], label = grp_label, ...) +
      guides(color = guide_legend(title=cell_label,override.aes = list(stroke = 1, 
                                                                       alpha = 1, shape = 16, size = pt_size, color=grp_color[1:ngrp])), 
             alpha =  "none")
  }
  
  require(Seurat)
  seu <- CreateSeuratObject(counts=corMat)
  Idents(seu) <- factor(cluster_orderd, levels=1: max(cluster_orderd))
  
  seu[["RNA"]]@scale.data <- corMat
  
  doHeatmap(seu, features = row.names(seu), cell_label=legend_title, pt_size=pt_size,
            grp_color=grp_color, grp_label=grp_label)
  
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


theme_hip2 <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                        base_rect_size = base_size/22, border_color = 'gray10', fill_color='gray80')
{
  # panel.background = element_rect(fill = fill_color,
  #                                 colour = fill_color),
  half_line <- base_size/2
  t <- theme(panel.border = element_rect(fill = NA,
                                         colour = border_color),
             #line = element_blank(), #rect = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE,
             panel.spacing = unit(half_line/2, "pt"), plot.margin = unit(rep(0.2,4), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}

theme_bulb16 <- function (base_size = 11, base_family = "", base_line_size = base_size/22,
                        base_rect_size = base_size/22, border_color = 'gray10', fill_color='gray80')
{
  # panel.background = element_rect(fill = fill_color,
  #                                 colour = fill_color),
  half_line <- base_size/2
  t <- theme(panel.border = element_rect(fill = NA,
                                         colour = border_color),
             #line = element_blank(), #rect = element_blank(),
             text = element_text(family = base_family, face = "plain",
                                 colour = "black", size = base_size, lineheight = 0.9,
                                 hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
                                 debug = FALSE), axis.text = element_blank(), axis.title = element_blank(),
             axis.ticks.length = unit(0, "pt"), axis.ticks.length.x = NULL,
             axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL,
             axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL,
             axis.ticks.length.y.right = NULL, legend.box = NULL,
             legend.key.size = unit(1.2, "lines"), legend.position = "right",
             legend.text = element_text(size = rel(0.8)), legend.title = element_text(hjust = 0),
             strip.text = element_text(size = rel(0.8)), strip.switch.pad.grid = unit(half_line/2,
                                                                                      "pt"), strip.switch.pad.wrap = unit(half_line/2,
                                                                                                                          "pt"), panel.ontop = FALSE,
             panel.spacing = unit(half_line/2, "pt"), plot.margin = unit(rep(0.2,4), "lines"),
             plot.title = element_text(size = rel(1.2), hjust = 0,
                                       vjust = 1, margin = margin(t = half_line)), plot.title.position = "panel",
             plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(t = half_line)),
             plot.caption = element_text(size = rel(0.8), hjust = 1,
                                         vjust = 1, margin = margin(t = half_line)), plot.caption.position = "panel",
             plot.tag = element_text(size = rel(1.2), hjust = 0.5,
                                     vjust = 0.5), plot.tag.position = "topleft",
             complete = TRUE)
  #ggplot2:::ggplot_global$theme_all_null %+replace% t
  t + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}

myRidgePlot <- function (object, features,group.names='Domain' , cols = NULL, idents = NULL, sort = FALSE, 
                         assay = NULL, group.by = NULL, y.max = NULL, same.y.lims = FALSE, 
                         log = FALSE, ncol = NULL, slot = "data", stack = FALSE, 
                         combine = TRUE, fill.by = "feature") {
  require(ggplot2)
  return(ExIPlot(object = object,group.names =group.names, type = "ridge", features = features, 
                 idents = idents, ncol = ncol, sort = sort, assay = assay, 
                 y.max = y.max, same.y.lims = same.y.lims, cols = cols, 
                 group.by = group.by, log = log, slot = slot, stack = stack, 
                 combine = combine, fill.by = fill.by))
}
#ExIPlot(seu2, features = DEGList[[k]][1], type='ridge') + ylab("Domain")
ExIPlot <- function (object, features, group.names='Domain', type = "violin", idents = NULL, 
          ncol = NULL, sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE, 
          adjust = 1, cols = NULL, pt.size = 0, group.by = NULL, split.by = NULL, 
          log = FALSE, slot = "data", stack = FALSE, combine = TRUE, 
          fill.by = NULL, flip = FALSE, raster = NULL, title_size=20) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning("'ncol' is ignored with 'stack' is TRUE", 
              call. = FALSE, immediate. = TRUE)
    }
    if (!is.null(x = y.max)) {
      warning("'y.max' is ignored when 'stack' is TRUE", 
              call. = FALSE, immediate. = TRUE)
    }
  }
  else {
    ncol <- ncol %||% ifelse(test = length(x = features) > 
                               9, yes = 4, no = min(length(x = features), 3))
  }
  data <- FetchData(object = object, vars = features, slot = slot)
  pt.size <- pt.size %||% AutoPointSize(data = object)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  }
  else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% 
                                                 idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  }
  else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  }
  else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    else if (length(x = cols) == 1 && cols == "interaction") {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    }
    else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiExIPlot(type = type, data = data, idents = idents, 
                        split = split, sort = sort, same.y.lims = same.y.lims, 
                        adjust = adjust, cols = cols, pt.size = pt.size, 
                        log = log, fill.by = fill.by, flip = flip))
  }
  if(type =='ridge'){
    plots <- lapply(X = features, FUN = function(x) {
      return(SingleExIPlot(type = type, data = data[, x, drop = FALSE], 
                           idents = idents, split = split, sort = sort, y.max = y.max, 
                           adjust = adjust, cols = cols, pt.size = pt.size, 
                           log = log, raster = raster)+ ggplot2::ylab(group.names))
    })
  }else{
    plots <- lapply(X = features, FUN = function(x) {
      return(SingleExIPlot(type = type, data = data[, x, drop = FALSE], 
                           idents = idents, split = split, sort = sort, y.max = y.max, 
                           adjust = adjust, cols = cols, pt.size = pt.size, 
                           log = log, raster = raster) + ggplot2::ylab(group.names))
    })
  }
    
  label.fxn <- switch(EXPR = type, violin = if (stack) {
    xlab
  } else {
    ylab
  }, splitViolin = if (stack) {
    xlab
  } else {
    ylab
  }, ridge = xlab, stop("Unknown ExIPlot type ", type, 
                        call. = FALSE))
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1], 
                  "_")
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = "DimReduc")) {
        plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
      }
      else if (inherits(x = object[[obj]], what = "Assay")) {
        next
      }
      else {
        warning("Unknown object type ", class(x = object), 
                immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  for (i in 1:length(x = plots)) {
    plots[[i]] <- plots[[i]] +  ggplot2::ylab(group.names)
  }
  
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

featurePlot <- function(seu, feature, cols, pt_size, quant = 0.9){
  dat <- as.data.frame(seu[["Spatial"]]@cell.embeddings)
  dat$Expression <- seu[['RNA']]@scale.data[toupper(feature),]
  #med <- median(seu[['RNA']]@scale.data[toupper(feature),])
  med <- quantile(seu[['RNA']]@scale.data[toupper(feature),], quant)
  ggplot(data=dat, aes(x=Spatial_1, y=Spatial_2, color=Expression)) + geom_point(size=pt_size) +
    #scale_colour_gradient(low = cols[1],high = cols[2]) + 
    scale_colour_gradient2(low = "#0571B099", mid = "white", high = "#CA0020", midpoint = med)+
    #scale_colour_gradient(low = "#0571B0",high = "#CA0020") + 
    mytheme_graybox() + 
    ggtitle(feature) + theme(title =element_text(size=title_size, color=1, face='italic'), legend.position = 'none')
  
}

featurePlot_pseudoassoc <- function(seu, feature,  pt_size, quant = 0.9){
  dat <- as.data.frame(seu[["Spatial"]]@cell.embeddings)
  dat$Expression <- seu[['RNA']]@scale.data[toupper(feature),]
  #med <- median(seu[['RNA']]@scale.data[toupper(feature),])
  med <- quantile(seu[['RNA']]@scale.data[toupper(feature),], probs =quant)
  ggplot(data=dat, aes(x=Spatial_1, y=Spatial_2, color=Expression)) + geom_point(size=pt_size) +
    #scale_colour_gradient(low = cols[1],high = cols[2]) + 
    #scale_colour_gradient2(low = "#0571B099", mid = "white", high = "#CA0020", midpoint = med)+
    #scale_colour_gradient(low = "#0571B0",high = "#CA0020") + 
    scale_colour_gradient2(low = "#26986D", mid = "#F8F7DC", high = "#CA0020", midpoint = med)+
    mytheme_graybox() + 
    ggtitle(feature) + theme(title =element_text(size=title_size, color=1, face='italic'), legend.position = 'none')
}



scaleFUN <- function(x) sprintf("%.2f", x)

bar_plot <- function(percentage_long,geom_bar_position,legend_position,color_pal){
  ggplot(percentage_long, aes(y = value, x = factor(Cluster), fill = variable)) +        ## global aes
    scale_fill_manual(values= color_pal,name = 'Cell Type')+
    scale_y_continuous(labels=scaleFUN) +
    geom_bar(position=geom_bar_position, stat="identity",width=0.7,color="black") +
    ggtitle(paste("")) +
    theme_bw()+xlab("")+ylab("")+
    theme_classic() +
    theme(plot.title = element_text(size = 14,hjust = 0.5),
          text = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.line = element_line(colour = "grey"),
          legend.position = legend_position,
          panel.background = element_rect(fill = "white", colour = NA), panel.border = element_rect(fill = NA, colour = "grey"))
}
