## --- Utilities ------------------------------------------------------------ ##

# Introduction
# Contains utility and other custom functions

## --- Source other package utilities --------------------------------------- ##

## --- Load libraries ------------------------------------------------------- ##

library(tidyr)
library(dplyr)
library(homologene)
library(stringr)
library(gprofiler2)
library(SCpubr)
library(ggplot2)
library(Seurat)
library(cowplot)
library(SPOTlight)
library(ggcorrplot)
library(reshape2)
library(ggpubr)
library(viridis)
library(scatterpie)
library(ggrepel)
library(tibble)
library(stringdist)
library(ggbreak)
library(ggvenn)
library(harmony)
library(enrichR)
library(igraph)
library(visNetwork)
library(heatmaply)
library(scales)
library(ggsci)

# --- createEmptyObject -------------------------------------------------- ##

createEmptyObject <- function(object,
                              name = "object",
                              assay = "RNA",
                              subset = F,
                              subset_n = 60000) {

  #' Create an emtpy seurat object
  #'
  #' createEmptyObject, creates a new seurat object from an existing one, with
  #' the metadata and raw counts but minus all other downstream analysis
  #'
  #' @param object Seurat Object
  #' @param name Character String. Custom name of object.
  #' @param subset Logical. Whether to subset down the number of cells
  #' @param subset_n Numeric. What number N to subset
  #'
  #' @importFrom Seurat CreateSeuratObject
  #'
  #' @export


  # Create new seurat object with new genes and project.name and old meta.data
  new_object <- CreateSeuratObject(data, project = name, meta.data = meta)

  # Create new assay to differentiate normalized and unnormalized counts
  new_object@assays$LOG <- new_object@assays$RNA

  rm(data)

  if (save == T){

    dir2 <- paste(paste(save_dir, name, sep=""), ".rds", sep = "")

    if (sample == T) {

      dir2 <- paste(paste(save_dir, name, sep=""), "_sample_", N, ".rds", sep = "")

    }

    saveRDS(new_object, file = dir2)

  }

  if (return == T) {

    return(new_object)
  }


}

## --- doSpatialFeaturPlot -------------------------------------------------- ##

doSpatialFeaturePlot <- function(X,
                                 feature,
                                 assay = 'SCT',
                                 slot = 'data',
                                 subset = F,
                                 subset_meta = 'region',
                                 subset_ident = 'endocrine',
                                 min_val = 'default',
                                 max_val = 'default',
                                 scale = 'viridis',
                                 cols = NA) {


  # --- Documentation ------------------------------------------------------- ##

  #' @title doSpatialFeaturePlot
  #'
  #' @description This is a function to do a SpatialFeature plot
  #'
  #' @param X Seurat object
  #' @param features Vector of character strings. Genes and/or meta data to plot
  #' @param assay Character string. Name of assay to pull genes from. Default is "SCT"
  #' @param subset Logical. Whether to subset the object before plotting, set by subset_meta and subset_ident. Default is F
  #' @param slot Character string. Name of slot to pull genes from. Default is "data"
  #' @param default_col Logical. Whether to plot with default colors, i.e. rainbow. Default is T,
  #' @param cols Vector of character strings. Colors for plot in order of features to plot.
  #' @param filter_val Numeric. Values in below cutoff is treated as 0s.
  #' @param subset_meta Character string. Name of metadata to subset if subset == T. Default is 'region'
  #' @param subset_ident Character string. Name of ident to subset from subset_meta if subset == T. Default is 'endocrine'
  #' @param min_val,max_val Float. Limits for plot, options include "default" for automatic ggplot2 limits, "object_val" for the max value of the object, incase a subset is being used, or entere a numeric value. Default is "default"
  #' @param cols Vector of colors. Only for discrete scale.
  #'
  #' @export
  #'
  #' @import dplyr
  #' @import Seurat
  #' @import ggplot2
  #' @importFrom scatterpie geom_scatterpie

  # Get meta data in features

  object <- X

  meta_object <- FetchData(object, vars = c(feature))

  colnames(meta_object) <- "feature"

  if (subset == T) {

    object <- SetIdent(object, value = subset_meta)

    object <- subset(object, idents = subset_ident)

  }

  meta_subset <- FetchData(object, vars = c(feature))

  colnames(meta_subset) <- "feature"

  meta_subset$row <- object@images$slice1@coordinates$imagerow

  meta_subset$col <- object@images$slice1@coordinates$imagecol

  if (min_val == 'default') {

    min_val <- min(meta_subset$feature)

  } else if (min_val == 'object_val') {

    min_val <- min(meta_object$feature)

  }

  if (max_val == 'default') {

    max_val <- max(meta_subset$feature)

  } else if (max_val == 'object_val') {

    max_val <- max(meta_object$feature)

  } else {

    max_val <- max_val

    min_val <- min_val
  }


  if (scale == 'viridis') {

    p <- ggplot(meta_subset) +
      geom_point(aes(x=row,
                     y=col,
                     col = feature)) +
      scale_color_viridis_c(option = 'G',
                            limits = c(min_val, max_val), na.value = 'gray12') +
      theme_void() +
      theme(
        panel.border = element_rect(color = "gray12", fill = NA),
        panel.background = element_rect(fill = "gray12", color = "gray12"),
        plot.background = element_rect(fill = "gray12", color = "gray12"),
        legend.background = element_rect(fill = "gray12", color = "gray12"),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text(size=8, face = "bold", colour = "white", angle = 0),
        legend.text = element_text(size=8, face = "bold", colour = "white"),
        legend.title.align = 0.5) +
      coord_fixed()


  } else if (scale == 'gradient') {
    
    colorbar_breaks <- c(min_val, 0, max_val)
    colorbar_values <- rescale(colorbar_breaks)

    p <- ggplot(meta_subset) +
      geom_point(aes(x=row,
                     y=col,
                     col = feature)) +
      scale_color_gradient2(low = "#0000FF",
                            mid = 'gray12',
                            high = "#FF0000", midpoint = 0,
                            limits = c(min_val, max_val), na.value = 'gray12') +
      theme_void() +
      theme(
        panel.border = element_rect(color = "gray12", fill = NA),
        panel.background = element_rect(fill = "gray12", color = "gray12"),
        plot.background = element_rect(fill = "gray12", color = "gray12"),
        legend.background = element_rect(fill = "gray12"),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text(size=15, face = "bold", colour = "white", angle = 0),
        legend.text = element_text(size=15, face = "bold", colour = "white"),
        legend.title.align = 0.5) +
      coord_fixed()

  } else if (scale == 'brewer') {

    p <- ggplot(meta_subset) +
      geom_point(aes(x=row,
                     y=col,
                     col = feature)) +
      scale_color_distiller(direction = -1, palette = 'Greys', limits = c(min_val, max_val), na.value = 'black') +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "black"),
        plot.background = element_rect(fill = "black"),
        legend.background = element_rect(fill = "black"),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text(size=15, face = "bold", colour = "white", angle = 0),
        legend.text = element_text(size=15, face = "bold", colour = "white"),
        legend.title.align = 0.5) +
      coord_fixed()


  } else if (scale =='manual') {

    p <- ggplot(meta_subset) +
      geom_point(aes(x=row,
                     y=col,
                     col = feature)) +
      scale_color_manual(values = cols) +
      theme_void() +
      theme(
        panel.border = element_rect(color = "gray12", fill = NA),
        panel.background = element_rect(fill = "gray12"),
        plot.background = element_rect(fill = "gray12"),
        legend.background = element_rect(fill = "gray12"),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text(size=15, face = "bold", colour = "white", angle = 0),
        legend.text = element_text(size=15, face = "bold", colour = "white"),
        legend.title.align = 0.5) +
      coord_fixed()
  }


  return(p)

}

## --- doScatterPie --------------------------------------------------------- ##

# Adapted from SPOTlight https://github.com/MarcElosua/SPOTlight

doScatterpie <- function(X,
                         features,
                         assay = "SCT",
                         slot = "data",
                         default_col = T,
                         cols = NA,
                         subset = F,
                         filter = F,
                         filter_val = 0,
                         subset_meta = 'region',
                         subset_ident = 'endocrine')
{

  # --- Documentation ------------------------------------------------------- ##

  #' @title doScatterPie
  #'
  #' @description This is a function to do a scatterpie plot
  #'
  #' @param X Seurat object
  #' @param features Vector of character strings. Genes and/or meta data to plot
  #' @param assay Character string. Name of assay to pull genes from. Default is "SCT"
  #' @param subset Logical. Whether to subset the object before plotting, set by subset_meta and subset_ident. Default is F
  #' @param filter Logical. Whether to filter out values below a threshold set by filter_val. Default is F
  #' @param slot Character string. Name of slot to pull genes from. Default is "data"
  #' @param default_col Logical. Whether to plot with default colors, i.e. rainbow. Default is T,
  #' @param cols Vector of character strings. Colors for plot in order of features to plot.
  #' @param filter_val Numeric. Values in below cutoff is treated as 0s.
  #' @param subset_meta Character string. Name of metadata to subset if subset == T. Default is 'region'
  #' @param subset_ident Character string. Name of ident to subset from subset_meta if subset == T. Default is 'endocrine'
  #'
  #' @export
  #'
  #' @import dplyr
  #' @import Seurat
  #' @import ggplot2
  #' @importFrom scatterpie geom_scatterpie

  # --- End of Documentation ------------------------------------------------ ##

  # Get meta data in features

  object <- X

  if (subset == T) {

    object <- SetIdent(object, value = subset_meta)

    object <- subset(object, idents = subset_ident)

  }

  meta <- FetchData(object, vars = c(features))

  if (filter == T) {

    meta[meta < filter_val] <- 0

  }

  meta$row <- object@images$slice1@coordinates$imagerow

  meta$col <- object@images$slice1@coordinates$imagecol


  features_final <- intersect(features, colnames(meta))

  if (default_col == T) {

    color <- rainbow(length(features_final))

  } else {

    color <- cols
  }

  names(color) <- features_final

  # Plot
  p  <- ggplot() + geom_scatterpie(
    data = meta,
    aes(
      x = row,
      y = col
    ),
    cols = features_final,
    color = NA,
    alpha = 1,
    pie_scale = 0.4) +
    scale_fill_manual(
      values = color,
      breaks = names(color)) +
    coord_fixed() +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      legend.background = element_rect(fill = "black"),
      axis.line=element_blank(),axis.text.x=element_blank(),
      axis.text.y=element_blank(),axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="bottom",
      legend.text = element_text(colour = 'white')
    )

  return(p)
}

## --- VlnDotPlot ----------------------------------------------------------- ##

VlnDotPlot <- function(object,
                       feature,
                       group,
                       fill_name = 'Average Score',
                       y_name = 'Absolute Score',
                       order = F) {

  data <- FetchData(object, vars = c(feature, group)) %>%
    rename_all(~ c('feature', 'group'))

  data_means <- data %>%
    group_by(group) %>%
    dplyr::summarise(mean_feature = mean(feature))

  data <- data %>%
    left_join(data_means, by = join_by(group == group))

  if (order == T) {

    group_order <- data_means %>%
      arrange(mean_feature) %>%
      pull(group)


  } else {

    group_order <- data_means %>%
      pull(group)

  }

  p <- ggplot(data, aes(x = factor(group, levels = group_order),
                        y = feature,
                        fill = mean_feature)) +
    geom_violin()  +
    scale_fill_viridis(na.value = 'black',
                       option = 'D',
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black",
                                              title = 'Avg. Score')) +
    stat_summary(aes(x = group,
                     y = feature), fun="mean",
                 geom="point", color="red", size = 1) +

    theme_classic() +
    theme( plot.title = element_blank(),
           axis.title.y = element_blank(),
           axis.title.x = element_text(size=8, face = "bold", colour = "black"),
           axis.text.x = element_text(size =8, angle = 90, face = "bold", colour = "black"),
           axis.text.y = element_text(size =8, face = "bold", colour = "black"),
           legend.title = element_text(size=8, face = "bold", colour = "black", angle = 0),
           legend.text = element_text(size=8, face = "bold", colour = "black"),
           legend.title.align = 0.5,
           legend.key.width = unit(0.5, 'cm'),
           legend.key.height = unit(0.3, 'cm'))+
    labs(fill = fill_name,
         y = y_name) +
    coord_flip()

  return(p)


}

## --- doCorrPlot --------------------------------------------------------- ##

doCorrPlot <- function(object,
                       features,
                       assay = "SCT",
                       slot = "data",
                       test = 'pearson',
                       cols = c("#6D9EC1", "white", "#E46726"),
                       rval = F,
                       order = T,
                       pval = F,
                       plottype = 'full') {

  # --- Documentation ------------------------------------------------------- ##

  #' @title doCorrPlot
  #'
  #' @description This is a function to do a correlation plot
  #'
  #' @param object Seurat object
  #' @param features Vector of character strings. Genes and/or meta data to plot
  #' @param assay Character string. Name of assay to pull genes from. Default is "SCT"
  #' @param slot Character string. Name of slot to pull genes from. Default is "data"
  #' @param cols Vector of character strings. Colors for correlation plot. Default is c("#6D9EC1", "white", "#E46726").
  #' @param rval Logical. Whether to plot pearsons R values in plot. Default is T
  #' @param order Logical. Whether to order the correlation matrix. Default is T.
  #' @param plottype Character string. Type of correlation matrix, full, lower or upper. Default is "full'.
  #' @param return Logical. Whether to return the plot object. Default is F
  #'
  #' @export
  #'
  #' @import dplyr
  #' @importFrom ggcorrplot ggcorrplot
  #' @importFrom stringr str_to_title

  # --- End of Documentation ------------------------------------------------ ##

  # Get meta data in features
  meta <- object@meta.data %>%
    select(any_of(features))

  # Get genes in features
  genes <- intersect(features, rownames(object))

  data <- GetAssayData(object,
                       assay = assay)[genes,] %>%
    as_tibble()

  # Transpose matrix because coercing dataframe with 1 column vs more is different
  if (ncol(data) == ncol(object)) {

    rownames(data) <- genes

    data <- t(data)

  } else {

    colnames(data) <- genes
  }


  # Merge gene and metadata information
  data_meta <- as_tibble(cbind(meta, data)) %>%
    as.matrix()

  # Get correlation
  corr <- cor(data_meta, method = test)

  if (pval == T) {

    # Get p-values
    p.mat <- cor_pmat(
      x = data_meta,
      conf_int = 0.95,
      method = test)

  } else {

    p.mat <- NULL

  }

  # Plot correlation matrix
  plot <- ggcorrplot(
    corr = corr,
    p.mat = p.mat,
    hc.order = order,
    lab = rval,
    colors = cols, type = plottype) +
    theme_classic() +
    theme(
      axis.line=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x = element_text(size = 8, angle = 60, vjust = 1, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 8, face = "bold"),
      legend.title = element_text(size=8, face = "bold"),
      legend.text = element_text(size=8, vjust = 1, hjust = 1, face = "bold"),
      legend.key.width = unit(0.3, 'cm'),
      legend.key.height = unit(0.2, 'cm')) +
    guides(fill = guide_colorbar(title = "Correlation"))


    return(plot)

}




## --- StackedVlnPlot -------------------------------------------------------- ##


## StackedVlnPlot Function
## https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
## (c) Ming Tang 2018
## Modified for visual effects

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1.5), angle = 0, face = "bold", colour = "black"),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {

  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle=90, face = "bold", colour = "black", size = rel(1.5)), axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


## --- filterMarkers -------------------------------------------------------- ##


filterMarkers <- function(markers,
                       gen_num = 10,
                       goi = "",
                       fold_cutoff = 0,
                       pval_cutoff = 0.05,
                       pct_cutoff = 0.1,
                       group_var = "cluster",
                       direction = 'both',
                       add_metrics = T,
                       make_unique = F,
                       unique_col = 'avg_log2FC') {

  #' Documentation
  #' @title filterMarkers
  #'
  #' @description This is a function that takes a FindMarkers/FindAllMarkers output and filters it based on certain criteria and adds additional metrics
  #'
  #' @param markers Dataframe. Output of FindMarkers and FindAllMarkers
  #' @param filter Logical. Whether to filter or not. Deafult is T
  #' @param gen_num Numeric or Character. Number of genes to output if "all", then the function returns all genes that meets cutoffs. Default is 10
  #' @param fold_cutoff Float. Value of log2FC cutoff, used as negative foldcutoff for negative genes. Default is 0
  #' @param pval_cutoff Float. Value of adjusted pvalue cutoff. Default is 0.25
  #' @param pct_cutff Float. Value of percentage expression in population 1 (pct.1) cutoff. Default is 0.1
  #' @param add_metrics Logical. Whether to add additional metrics, pct.diff and FDR_pvalue. Default is T
  #' @param group_var Character string. Grouping variable, if multiple groups included in markers. Defualt is 'cluster'
  #' @param direction Character string. Defines which direction of genes to select. "top" for Top/positive log2FC genes or "bottom" for bottom/negative log2FC genes or "both" for both top and bottom markers. Default is 'both'
  #'
  #' @import dplyr
  #' @import tidyr
  #'
  #' @export

  # Filter based on cutoffs

  filt_markers <- markers

  filt_markers <- markers %>%
    filter(p_val_adj < pval_cutoff) %>%
    filter(pct.1 > pct_cutoff)

  if ((group_var %in% colnames(markers)) == TRUE) {

    message("Grouping markers by ", group_var)

    if (direction == "both" ) {

      filt_markers_down <- filt_markers %>%
        group_by(cluster) %>%
        slice_min(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC < -fold_cutoff)

      filt_markers_up <- filt_markers %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC > fold_cutoff)

      filt_markers <- bind_rows(filt_markers_down, filt_markers_up) %>%
        arrange(desc(avg_log2FC))


    } else if (direction == "down") {

      filt_markers <- filt_markers %>%
        group_by(cluster) %>%
        slice_min(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC < -fold_cutoff) %>%
        arrange(avg_log2FC)

    } else if (direction == "up") {

      filt_markers <- filt_markers %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC > fold_cutoff) %>%
        arrange(desc(avg_log2FC))
    }

  } else {

    message("Grouping not found, no grouping of markers")

    if (direction == "both" ) {

      filt_markers_down <- filt_markers %>%
        slice_min(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC < -fold_cutoff)

      filt_markers_up <- filt_markers %>%
        slice_max(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC > fold_cutoff)

      filt_markers <- rbind(filt_markers_down, filt_markers_up) %>%
        arrange(avg_log2FC)


    } else if (direction == "bottom") {

      filt_markers <- filt_markers %>%
        slice_min(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC < -fold_cutoff) %>%
        arrange(desc(avg_log2FC))

    } else if (direction == "top") {

      filt_markers <- filt_markers %>%
        slice_max(order_by = avg_log2FC, n = gen_num) %>%
        filter(avg_log2FC > fold_cutoff) %>%
        arrange(avg_log2FC)
    }
  }


  if (add_metrics == T) {

    filt_markers <- filt_markers %>%
      mutate(
        # Percentage diff allows detection of lower expressed genes
        pct.diff = pct.1 - pct.2,
        # Weighted percentage diff allows detection of lower expressed genes
        pct.weight.diff = 1 - ((1 / pct.1) * pct.2),
        # Bonferroni log
        log_pval_bon = -log(p_val_adj),
        # FDR adjusted p_val
        p_val_adj_bh = p.adjust(p_val, method="BH"),
        # FDR log
        log_pval_bh = -log(p_val_adj_bh),
        # Rank log2FC
        rank_FC_pctdiff = rank(avg_log2FC*pct.diff),
        rank_FC_pctweight = rank(avg_log2FC*pct.weight.diff),
        rank_FC_pval_bon = rank(avg_log2FC*log_pval_bon),
        rank_FC_pval_bh = rank(avg_log2FC*log_pval_bh)
      )

  }

  if (('gene' %in% colnames(filt_markers)) == FALSE) {

    filt_markers$gene <- rownames(filt_markers)

  }

  # Filter based on input of genes
  if (length(goi) > 1) {

    goi_int <- which(filt_markers$gene %in% goi)

    filt_markers <- filt_markers[goi_int,]

  }

  if(make_unique == T) {

    filt_markers <- filt_markers %>%
      group_by(gene) %>%
      arrange(desc(unique_col)) %>%
      slice_head(n = 1)

  }
  return(filt_markers)

}

## --- convertGenes --------------------------------------------------------- ##
# --- convertGenes

# --- Documentation ----------------------------------------------------- ##

#' convertGenes
#'
#'This is a function to convert mouse to human genes or vice verse using
#'the homologene package
#'
#'@param genes Vector variable of genes to convert
#'@param converTo Character string variable specifying which genome to convert
#'to either "human", i.e convert mouse to human genes, or "mouse", i.e. convert 
#' human to mouse genes.
#'@param caseChange Logical variable, if True, genes not converted are uppercased
#'or lowercased based on genome to convert
#'@param match Logical variable, if True, converted genes are outputted to match 
#'indice of input genes variable. If true, it removes all duplicated genes that 
#'have different converted genes. Set to false if all genes should be outputted
#'

# --- End of Documentation ---------------------------------------------- ##


convertGenes_deprecated <- function(genes, 
                         convertTo = "human", 
                         caseChange = T, 
                         match = T) {
  
  require(homologene)
  require(stringr)
  
  if (convertTo == "human") {
    
    genesTable <- mouse2human(genes)
    
    genesTable <- genesTable[,1:2]
    
    colnames(genesTable) = c("Mouse", "Human")
    
    genesM <- genesTable[,1]
    
    genesH <- genesTable[,2]
    
    
    if (caseChange == T) {
      
      remGenes <- setdiff(genes, genesM)
      
      remGenes_upper <- toupper(remGenes)
      
      remGenes_df <- data.frame("Mouse" = remGenes, "Human" = remGenes_upper)
      
      genesTable <- rbind(genesTable, remGenes_df)
      
    }
    
    genesM <- genesTable[,1]
    
    genesH <- genesTable[,2]
    
    
    if (any(duplicated(genesM)) == T) {
      
      indDupM <- duplicated(genesM) |duplicated(genesM, fromLast = T)
      
      print("The following mouse genes have multiple human homologues")
      
      print(genesTable[indDupM,])
      
    }
    
    if (any(duplicated(genesH)) == T) {
      
      indDupH <- duplicated(genesH) |duplicated(genesH, fromLast = T)
      
      print("The following mouse genes have the same human homologues")
      
      print(genesTable[indDupH,])
      
    }
    
    if (match == T) {
      
      if (any(duplicated(genesM)) == T) {
        
        indDupM <- duplicated(genesM) 
        
        likelyConvert <- toupper(genesM[indDupM])
        
        presentGenes <- likelyConvert[na.omit(match(genesH,likelyConvert))]
        
        if (length(presentGenes) > 0 ) {
          
          indPresent <-  match(str_to_title(presentGenes),genesM)
          
          genesTable[indPresent, "Human"] <- presentGenes
          
        }
        
        genesH <- genesTable[,2]
        
        print("Multiple unique homologues of the following genes deleted except best fit")
        
        print(unique(genesTable[indDupM,"Mouse"]))
        
      }
      
      genesH <- genesH[match(genes, genesM)]
      
    }
    
    return (genesH)
    
  }  else if (convertTo == "mouse") {
    
    genesTable <- human2mouse(genes)
    
    genesTable <- genesTable[,1:2]
    
    colnames(genesTable) = c("Human", "Mouse")
    
    genesH <- genesTable[,1]
    
    genesM <- genesTable[,2]
    
    
    if (caseChange == T) {
      
      remGenes <- setdiff(genes, genesH)
      
      remGenes_upper <- str_to_title(remGenes)
      
      remGenes_df <- data.frame("Human" = remGenes, "Mouse" = remGenes_upper)
      
      genesTable <- rbind(genesTable, remGenes_df)
      
    }
    
    genesH <- genesTable[,1]
    
    genesM <- genesTable[,2]
    
    
    if (any(duplicated(genesH)) == T) {
      
      indDupH <- duplicated(genesH) |duplicated(genesH, fromLast = T)
      
      print("The following human genes have multiple mouse homologues")
      
      print(genesTable[indDupH,])
      
    }
    
    if (any(duplicated(genesM)) == T) {
      
      indDupM <- duplicated(genesM) |duplicated(genesM, fromLast = T)
      
      print("The following human genes have the same mouse homologues")
      
      print(genesTable[indDupM,])
      
    }
    
    if (match == T) {
      
      if (any(duplicated(genesH)) == T) {
        
        indDupH <- duplicated(genesH) 
        
        likelyConvert <- str_to_title(genesH[indDupH])
        
        presentGenes <- likelyConvert[na.omit(match(genesM,likelyConvert))]
        
        if (length(presentGenes) > 0 ) {
          
          indPresent <-  match(toupper(presentGenes),genesH)
          
          genesTable[indPresent, "Mouse"] <- presentGenes
          
        }
        
        genesM <- genesTable[,2]
        
        print("Multiple unique homologues of the following genes deleted except best fit")
        
        print(unique(genesTable[indDupH,"Human"]))
        
      }
      
      genesM <- genesM[match(genes, genesH)]
      
    }
    return (genesM)
  }
  
}

updateGenes_deprecated <- function(features,geneDB, truncate = T){
  
  new_features <- features
  
  if (truncate == T) {
    
    notupd_ind_feat <- which(!(features%in% geneDB[,1]))
    features_trunc <- features[notupd_ind_feat]
    
    
    notupd_ind_gdb <- which(!(geneDB[,1]%in% features))
    geneDB_trunc <- geneDB[notupd_ind_gdb,]
    
    features <- features_trunc
    
    geneDB <- geneDB_trunc
  }
  
  
  feature_list <- list()
  
  for (i in 1:nrow(geneDB)){
    
    vec <- as.vector(unlist(geneDB[i,]))
    
    vec_split <- unlist(strsplit(vec, ", "))
    
    in_ind <- features %in% vec
    
    if (any(in_ind) ==T) {
      
      gene <- features[in_ind] 
      
      rows <- names(feature_list)
      
      for (j in 1:length(gene)) {
        
        row_ind <- gene[j] %in% rows
        
        if (any(row_ind) == T){
          
          feature_list[[gene[j]]] <- c(feature_list[[gene[j]]], geneDB[i,1])
          
        } else {
          
          feature_list[gene[j]] <- geneDB[i,1]
          
        }
        
      }
    }
  } 
  
  lengths <- sapply(feature_list, length)
  
  sing_genes <- names(which(lengths==1))
  
  upd_features <- as.vector(unlist(feature_list[sing_genes]))
  
  sing_ind <- match(sing_genes, new_features)
  
  new_features[sing_ind] <- upd_features
  
  dup_genes <- names(which(lengths > 1))
  
  dups <- feature_list[dup_genes]
  
  dups <- lapply(dups, unique)
  
  if(length(dups) < 1) {
    
    cat("No genes have multiple aliases")
  } else {
    for (i in 1:length(dups)) {
      cat("These genes have multiple aliases", names(dups)[i], ": ")
      print(dups[[i]])
    }
  }
  return(list(new_features, dups))
  
} 

## --- convertGenes --------------------------------------------------------- ##

# convertGenes
convertGenes <- function(genes,
                         convert_to = "human") {

  # --- Documentation ------------------------------------------------------- ##

  #' @title convertGenes
  #'
  #' @description This is a function to convert mouse to human genes or vice verse using
  #' the homologene package
  #'
  #' @param genes Vector of character strings. Genes to convert
  #' @param conver_to Character string. Specifying which genome to convert
  #' to. Either of "human", i.e convert mouse to human genes, or "mouse", i.e. convert
  #' human to mouse genes.
  #' @param case_change Boolean. If True, genes not converted are uppercased
  #' or lowercased based on genome to convert
  #' @param remove_dups Boolean. If True, converted genes are outputted to match
  #' indice of input genes variable. If true, it removes all duplicated genes that
  #' have different converted genes. Set to false if all genes should be outputted
  #'
  #' @return A list of converted genes
  #'
  #' @export
  #'
  #' @importFrom homologene mouse2human
  #' @importFrom homologene human2mouse
  #' @importFrom stringr str_to_title

  # --- End of Documentation ------------------------------------------------ ##


  # Similarity function
  similarr <- function(x,y) {
    return(1 - stringdist(x, y, method = "lv")/ max(nchar(x), nchar(y)))
  }
  
  # Make new features variable
  features <- genes
  
  if (convert_to == "human") {
    
    genes_table <- mouse2human(features) # Convert using homologene package
    genes_table <- genes_table[,1:2] # Get only first two relevant columns
    colnames(genes_table) = c("original", "updated")
    
    # Remove ghost genes
    genes_table <- genes_table %>%
      filter(original %in% features)
    
    # Get absent features
    features_not_present <- setdiff(features, genes_table$original)
    
    # Make new dataframe with absent features
    unchaged_genes_df <- data.frame("original" = features_not_present,
                                    "updated" = toupper(features_not_present))
    
    # Bind to original genes_table
    genes_table <- rbind(genes_table, unchaged_genes_df)
    
    genes_table <- genes_table %>%
      mutate(cased = toupper(original), # Make a case changed column for comparisions
             alternate = cased, # Make an alternate column that will be switched based on iteration
             final = updated, # Make final column that will be assessed for duplicates
             duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F), # Get the first set of final duplicates
             duplicated_original = ifelse(duplicated(original) | duplicated(original, fromLast = T), T, F), # Get the first set of original duplicates
             similarity = similarr(final, cased)) # Get the similarity score between final and "original" or cased
    
  } else if (convert_to == "mouse") { # Repeat for mouse
    
    genes_table <- human2mouse(features) 
    genes_table <- genes_table[,1:2] 
    colnames(genes_table) = c("original", "updated")
    
    genes_table <- genes_table %>%
      filter(original %in% features)
    
    features_not_present <- setdiff(features, genes_table$original)
    
    unchaged_genes_df <- data.frame("original" = features_not_present,
                                    "updated" = str_to_title(features_not_present)) # str_to_title instead of toupper
    
    genes_table <- rbind(genes_table, unchaged_genes_df)
    
    genes_table <- genes_table %>%
      mutate(cased = str_to_title(original), # str_to_title instead of toupper
             alternate = cased,
             final = updated, 
             duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F),
             duplicated_original = ifelse(duplicated(original) | duplicated(original, fromLast = T), T, F),
             similarity = similarr(final, cased))
  }
  
  
  genes_table_orig_dup <- genes_table %>%
    filter(duplicated_original == F) # Get rid of original duplicates, which are handled later on
  
  finaldups <- any(genes_table_orig_dup$duplicated_final == T) # Check for final duplicates that are not original duplicates
  
  counter <- 1 # Counter for iteration
  
  while(finaldups && counter < 15) {
    
    # Get all duplicated final genes
    final_duplicates <- genes_table %>%
      filter(duplicated_final == T) %>%
      pull(final)
    
    # Split into two datasets
    # No duplicates in neither updated or cased (original)
    genes_table_no_dup <- genes_table %>%
      filter(!(updated %in% final_duplicates) & !(cased %in% final_duplicates))
    
    # With duplicates in either updated or cased (original)
    genes_table_dup <- genes_table %>%
      filter((updated %in% final_duplicates) | (cased %in% final_duplicates)) %>%
      mutate(similarity = similarr(final, cased)) %>% # recalculate similarity score with original
      group_by(final) %>% # group_by final
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # rank based on similarity score, window function works on groups
      ungroup() %>% # ungroup
      mutate(final = ifelse(rank == 1, final, alternate), # Keep final if similarity rank == 1, change to alternate if not
             alternate = ifelse(final == cased, updated, cased)) # Change alternate based on what final is
    
    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F)) # Get new duplicated finals
    
    
    # Check for duplicates in final again
    genes_table_orig_dup <- genes_table %>%
      filter(duplicated_original == F)
    
    finaldups <- any(genes_table_orig_dup$duplicated_final == T)
    
    counter <- counter + 1
    
  }
  
  if (any(genes_table$duplicated_original == T)) {
    # Now handle original duplicates
    # Get all original duplicate genes
    original_duplicates <- genes_table %>%
      filter(duplicated_original == T) %>%
      pull(original)
    
    # Split into two datasets based on duplicates
    genes_table_no_dup <- genes_table %>%
      filter(!(original %in% original_duplicates))
    
    genes_table_dup <- genes_table %>%
      filter(original %in% original_duplicates) %>%
      mutate(similarity = similarr(final, cased)) %>% # Get similarity score again 
      group_by(cased) %>% # Group by ORIGINAL/Cased
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # Rank similiarty
      slice_max(order_by = desc(rank)) %>% # Remove all other duplicates except first rank
      ungroup() # ungroup
    
    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F))
  }
  
  if (any(genes_table$duplicated_final == T)) {
    
    cat("Error duplicates still exist")
    
  }
  
  genes_table <- genes_table %>%
    arrange(match(original, features)) # Match the original order
  
  return(genes_table$final)
}

## --- updateGenes ---------------------------------------------------------- ##

updateGenes<- function(genes,
                       db,
                       print = F) {

  # --- Documentation ------------------------------------------------------- ##

  #' updateGenes
  #'
  #'This is a function to update all genes using gprofilers gconvert function
  #'
  #' @param genes Vector variable of genes to convert
  #' @param db Dataframe variable. HGNC symbol list.
  #' @param print Logical. Whether to print the updated genes or not
  #'
  #' @importFrom stringr str_split_fixed
  #' @importFrom stringdist stringdist
  #'
  #' @export

  # --- End of Documentation ------------------------------------------------ ##

  # Similarity function
  similarr <- function(x,y) {
    return(1 - stringdist(x, y, method = "lv")/ max(nchar(x), nchar(y)))
  }
  
  db <- hgnc
  
  # Make geneDB original variable
  geneDB_orig <- db
  
  # Get previous symbols
  previous_symbols <- str_split_fixed(geneDB_orig$Previous.symbols, pattern = ", ", n = Inf)
  
  # Get alias symbols
  alias_symbols <- str_split_fixed(geneDB_orig$Alias.symbols, pattern = ", ", n = Inf)
  
  # Merge symbols
  all_symbols <- cbind(previous_symbols, alias_symbols)
  
  # Change symbol names
  colnames(all_symbols) <- mapply(paste0, rep("original_symbol_",
                                              ncol(all_symbols)), c(1:ncol(all_symbols)))
  
  # Pivot longer
  geneDB <- geneDB_orig %>%
    select(Approved.symbol) %>% # Select only the approved column
    rename(updated = Approved.symbol) %>% # rename to "approved"
    cbind(all_symbols) %>% # cbind all the rest of the symbols
    mutate_all(~ ifelse(. == "", NA, .)) %>% # Make all "" into NAs
    pivot_longer(cols = starts_with("original"),
                 values_to = "original",
                 values_drop_na = T) %>% #pivot longer
    select(updated, original) # Select only approved and alternate columns
  
  
  # Make a new features variable
  features <- genes
  
  # Get approved features
  features_approved <- unique(intersect(features, geneDB$updated))
  
  # Get all genes that are not approved
  features_not_approved <- unique(setdiff(features, features_approved))
  
  # Get all genes that are not present
  features_not_present <- unique(setdiff(features_not_approved, geneDB$original))
  
  # Get all features that can't be changed
  features_unchaged <- unique(c(features_approved, features_not_present))
  
  # Get all genes that can be updated
  features_need_changing <- unique(setdiff(features_not_approved, features_not_present))
  
  # Remove ghost genes
  genes_table <- geneDB %>%
    filter(original %in% features_need_changing)
  
  # Make new dataframe with unchangeable features
  unchaged_genes_df <- data.frame("original" = features_unchaged,
                                  "updated" = features_unchaged) 

  
  # Bind to original genes_table
  genes_table <- rbind(genes_table, unchaged_genes_df)
  
  genes_table <- genes_table %>%
    mutate(alternate = original, # Make an alternate column that will be switched based on iteration
           final = updated, # Make final column that will be assessed for duplicates
           duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F), # Get the first set of final duplicates
           duplicated_original = ifelse(duplicated(original) | duplicated(original, fromLast = T), T, F), # Get the first set of original duplicates
           similarity = similarr(final, original)) # Get the similarity score between final and "original" or cased
  
  
  genes_table_orig_dup <- genes_table %>%
    filter(duplicated_original == F) # Get rid of original duplicates, which are handled later on
  
  finaldups <- any(genes_table_orig_dup$duplicated_final == T) # Check for final duplicates that are not original duplicates
  
  counter <- 1 # Counter for iteration
  
  while(finaldups && counter < 15) {
    
    # Get all duplicated final genes
    final_duplicates <- genes_table %>%
      filter(duplicated_final == T) %>%
      pull(final)
    
    # Split into two datasets
    # No final duplicates in neither updated or  original
    genes_table_no_dup <- genes_table %>%
      filter(!(updated %in% final_duplicates) & !(original %in% final_duplicates))
    
    # With final duplicates in either updated or original
    genes_table_dup <- genes_table %>%
      filter((updated %in% final_duplicates) | (original %in% final_duplicates)) %>%
      mutate(similarity = similarr(final, original)) %>% # recalculate similarity score with original
      group_by(final) %>% # group_by final
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # rank based on similarity score, window function works on groups
      ungroup() %>% # ungroup
      mutate(final = ifelse(rank == 1, final, alternate), # Keep final if similarity rank == 1, change to alternate if not
             alternate = ifelse(final == original, updated, original)) # Change alternate based on what final is
    
    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F)) # Get new duplicated finals
    
    
    # Check for duplicates in final again
    genes_table_orig_dup <- genes_table %>%
      filter(duplicated_original == F)
    
    finaldups <- any(genes_table_orig_dup$duplicated_final == T)
    
    counter <- counter + 1
    
  }
  
  if (any(genes_table$duplicated_original == T)) {
    
    # Now handle original duplicates
    # Get all original duplicate genes
    original_duplicates <- genes_table %>%
      filter(duplicated_original == T) %>%
      pull(original)
    
    # Split into two datasets based on duplicates
    genes_table_no_dup <- genes_table %>%
      filter(!(original %in% original_duplicates))
    
    genes_table_dup <- genes_table %>%
      filter(original %in% original_duplicates) %>%
      mutate(similarity = similarr(final, original)) %>% # Get similarity score again 
      group_by(original) %>% # Group by ORIGINAL/Cased
      mutate(rank = rank(desc(similarity), ties.method = "random")) %>%  # Rank similiarty
      slice_max(order_by = desc(rank)) %>% # Remove all other duplicates except first rank
      ungroup() # ungroup
    
    genes_table <- bind_rows(genes_table_no_dup, genes_table_dup) %>% # Rebind datasets
      mutate(duplicated_final = ifelse(duplicated(final) | duplicated(final, fromLast = T), T, F))
    
  }
  
  if (any(genes_table$duplicated_final == T)) {
    
    cat("Error duplicates still exist")
    
  }
  
  genes_table <- genes_table %>%
    arrange(match(original, features)) # Match the original order
  
  return(genes_table$final)
}





## --- removeDuplicates------------------------------------------------------ ##

removeDuplicates <- function(genes_table,
                             convert = 'upper') {

  # --- Documentation ------------------------------------------------------- ##

  #' removeDuplicates
  #'
  #' This is a function to remove duplicates in a gene_table of input and targets.
  #'
  #' @param genes_table Dataframe. Containing character strings of input and target gene
  #' @param convert Character string. Specifies whether to convert the first column to 'lower', lowercase using str_to_title or 'upper', uppercase using toupper from stringr. Default is upper
  #'
  #' @importFrom gprofiler2 gconvert
  #' @importFrom stringr toupper
  #' @importFrom stringr str_to_title
  #'
  #' @export

  # --- End of Documentation ------------------------------------------------ ##

  # Get input genes
  input_genes <- genes_table[,1]

  original_inputs <- input_genes

  # Get target genes
  target_genes <- genes_table[,2]

  # Make new dataframe of input and targets
  genes_table <- data.frame("Input" = input_genes,
                            "Target" = target_genes)

  # find all duplidcated input indices
  dup_input_ind <- sort(unique(c(which(duplicated(input_genes)),
                                         which(duplicated(input_genes, fromLast = T)))))

  # Find all duplicated input genes
  dup_input_genes <- unique(input_genes[dup_input_ind])


  # Find all duplicated target indices
  dup_target_ind <- sort(unique(c(which(duplicated(target_genes)),
                                  which(duplicated(target_genes, fromLast = T)))))

  # Find all duplicated target genes
  dup_target_genes <- unique(target_genes[dup_target_ind])


  # Get all original inputs that will need to be modifieid

  updated_inputs <- c()

  if (any(duplicated(input_genes))) {

    cat("The following input genes have multiple target genes\n")

    for (i in 1:length(dup_input_ind)) {

      cat(paste(genes_table[dup_input_ind[i], 1], " : ", genes_table[dup_input_ind[i], 2], sep = ""), "\n")

    }

  } else if (!any(duplicated(input_genes))) {

    cat("No input genes have multiple targets\n")

  }


  if (any(duplicated(target_genes))) {

    cat("\nThe following target genes have multiple input genes\n")

    for (i in 1:length(dup_target_ind)) {

      cat(paste(genes_table[dup_target_ind[i], 1], " : ", genes_table[dup_target_ind[i], 2], sep = ""), "\n")

    }


  } else if (!any(duplicated(target_genes))) {

    cat("\nNo target genes have multiple inputs\n")

  }

  if (any(duplicated(genes_table$Input))) {

    input_genes <- genes_table$Input

    target_genes <- genes_table$Target

    # find all duplidcated input indices
    dup_input_ind <- sort(unique(c(which(duplicated(input_genes)),
                                   which(duplicated(input_genes, fromLast = T)))))

    # Find all duplicated input genes
    dup_input_genes <- unique(input_genes[dup_input_ind])


    # Find all duplicated target indices
    dup_target_ind <- sort(unique(c(which(duplicated(target_genes)),
                                    which(duplicated(target_genes, fromLast = T)))))

    # Find all duplicated target genes
    dup_target_genes <- unique(target_genes[dup_target_ind])


    for (i in 1:length(dup_input_genes)) {

      gene <- dup_input_genes[i]

      gene_ind <- which(genes_table$Input == gene)

      # Get targets
      targets <- genes_table$Target[gene_ind]

      match_ind <- match(gene, targets)

      if (!is.na(match_ind)) {

        filt_ind <- gene_ind[-match_ind]

        genes_table <- genes_table[-filt_ind,]

        updated_inputs <- c(updated_inputs, gene)

      } else if (is.na(match_ind)) {

        genes_table <- genes_table[-gene_ind[-1],]

        updated_inputs <- c(updated_inputs, gene)

      }

    }

  }

  while (any(duplicated(genes_table$Target))) {

    target_genes <- genes_table$Target

    input_genes <- genes_table$Input

    dup_target_ind <- sort(unique(c(which(duplicated(target_genes)),
                                        which(duplicated(target_genes, fromLast = T)))))


    dup_target_genes <- unique(target_genes[dup_target_ind])


    for (i in 1:length(dup_target_genes)) {

      gene <- dup_target_genes[i]

      gene_ind <- which(genes_table$Target == gene)

      # Get inputs
      inputs <- genes_table$Input[gene_ind]

      match_ind <- match(gene, inputs)

      if (!is.na(match_ind)) {

        replace_inputs <- genes_table$Input[gene_ind[-match_ind]]

        genes_table$Target[gene_ind[-match_ind]] <- replace_inputs

        updated_inputs <- c(updated_inputs, replace_inputs)

      } else if (is.na(match_ind)) {

        replace_inputs <- genes_table$Input[gene_ind[-1]]

        genes_table$Target[gene_ind[-1]] <- replace_inputs

        updated_inputs <- c(updated_inputs, replace_inputs)

      }

    }

  }

  updated_ind <- unique(match(c(updated_inputs), genes_table$Input))

  cat("\nThe following conversions were made to handle duplicates\n")

  for (i in 1:length(updated_ind)) {

    cat(paste(genes_table[updated_ind[i], 1], " : ", genes_table[updated_ind[i], 2], sep = ""), "\n")

  }


  return(genes_table)

}

## --- RunVarimaxPCA -------------------------------------------------------- ##

RunVarimaxPCA <- function(object,
                          features = 'all',
                          reduction.name = "pca",
                          ncomp = 50,
                          norm = 'LOG') {

  #' Runs varimax rotation on PC
  #'
  #' Conducts varimax on the PCA loadings and then adds the new loadings to a different object
  #'
  #'@param object Seurat object with PCA already done
  #'@param reduction.name The name of the reduction the pca results is stored under
  #'@param ncomp Number of principle components to conduct varimax on
  #'@param norm Which assay to get scaled data from
  #'
  #'@importFrom stats varimax
  #'
  #'@return is a seurat object with varimax pca in a new reduction
  #'
  #'@export
  
  features <- intersect(features, rownames(object))


  pca_cell_embeddings <- object@reductions[[reduction.name]]@cell.embeddings
  pca_feature_loadings  <- object@reductions[[reduction.name]]@feature.loadings
  sdev <- object@reductions[[reduction.name]]@stdev
  
  if(any(features == 'all')) {
    data <- GetAssayData(object, assay = norm, slot = 'scale.data')
  } else {
    data <- FetchData(object, vars = features, slot = 'data')
  }


  data <- as.matrix(t(scale(data)))

  # prcomp_results$x is stored in cell.embeddings
  # prcomp_results$rotation is stored under feature.loadings
  # Modified from amoeba, 2013, stackoverflow

  rotatedLoadings <- varimax(pca_feature_loadings)$loadings
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  scores          <-  data %*% invLoadings

  #Convert to matrix (adapted from stackoverflow jay.sf, Dec, 2018)
  rotatedLoadings_matrix <- as.matrix(data.frame(matrix(as.numeric(rotatedLoadings),attributes(rotatedLoadings)$dim,dimnames=attributes(rotatedLoadings)$dimnames)))

  vpca <- CreateDimReducObject(
    embeddings = scores,
    loadings = rotatedLoadings_matrix,
    assay = assay,
    stdev = sdev,
    key = 'VPC_',
  )

  object@reductions[["vpca"]] <- vpca
  
  return(object)

}




## --- mergeMarkerDFs ------------------------------------------------------- ##

mergeMarkerDFs <- function(marker_list) {

  #' Documentation
  #' @title mergeMarkerDFs
  #'
  #' @description This is a function that merges a named list of marker dfs from the FindMarkres function, with the same column names
  #'
  #' @param marker_list List of named data frame. Output of FindMarkers and FindAllMarkers
  #'
  #' @import dplyr
  #' @import tidyr
  #'
  #' @export

  colnames <- setdiff(colnames(marker_list[[1]]), "updated_genes")


  for (i in 1:length(marker_list)) {

    marker_list[[i]]$dataset <- rep(names(marker_list[i]), nrow(marker_list[[i]]))

  }


  result <- bind_rows(marker_list) %>%
    pivot_wider(names_from = dataset, values_from = colnames, names_sep = "_")

}

## --- do_VolcanoPlot_modified ---------------------------------------------- ##

#' Compute a Volcano plot out of DE genes.
#'
#' @inheritParams doc_function
#' @param de_genes \strong{\code{\link[tibble]{tibble}}} | Output of `Seurat::FindMarkers()`.
#' @param pval_cutoff \strong{\code{\link[base]{numeric}}} | Cutoff for the p-value.
#' @param FC_cutoff \strong{\code{\link[base]{numeric}}} | Cutoff for the avg_log2FC.
#' @param plot_lines \strong{\code{\link[base]{logical}}} | Whether to plot the division lines.
#' @param line_color \strong{\code{\link[base]{character}}} | Color for the lines.
#' @param line_size \strong{\code{\link[base]{numeric}}} | Size of the lines in the plot.
#' @param add_gene_tags \strong{\code{\link[base]{logical}}} | Whether to plot the top genes.
#' @param add_tags_by \strong{\code{\link[base]{character}}} | Either "both", "positive" or "negative". To plot top genes on the negative side, positive side or both.
#' @param order_tags_by \strong{\code{\link[base]{character}}} | Either "both", "pvalue" or "logfc".
#' @param n_genes \strong{\code{\link[base]{numeric}}} | Number of top genes in each side to plot.
#' @param use_labels \strong{\code{\link[base]{logical}}} | Whether to use labels instead of text for the tags.
#' @param colors.use \strong{\code{\link[base]{character}}} | Color to generate a tetradic color scale with.
#'
#' @return A volcano plot as a ggplot2 object.
#' @export
#'
#' @example /man/examples/examples_do_VolcanoPlot.R
do_VolcanoPlot_modified <- function(sample,
                           de_genes,
                           pval_cutoff = 0.05,
                           FC_cutoff = 2,
                           pt.size = 2,
                           border.size = 1.5,
                           border.color = "black",
                           font.size = 14,
                           font.type = "sans",
                           plot.title = NULL,
                           plot.subtitle = NULL,
                           plot.caption = NULL,
                           plot_lines = TRUE,
                           line_color = "grey75",
                           line_size = 0.5,
                           add_gene_tags = TRUE,
                           add_tags_for = "both",
                           order_tags_by = "both",
                           n_genes = 5,
                           use_labels = FALSE,
                           colors.use = "steelblue",
                           plot.title.face = "bold",
                           plot.subtitle.face = "plain",
                           plot.caption.face = "italic",
                           axis.title.face = "bold",
                           axis.text.face = "plain",
                           legend.title.face = "bold",
                           legend.text.face = "plain"){
  # Add lengthy error messages.
  withr::local_options(.new = list("warning.length" = 8170))

  check_suggests(function_name = "do_VolcanoPlot")
  # Check if the sample provided is a Seurat object.
  check_Seurat(sample = sample)

  # Check logical parameters.
  logical_list <- list("add_gene_tags" = add_gene_tags,
                       "plot_lines" = plot_lines,
                       "use_labels" = use_labels)
  check_type(parameters = logical_list, required_type = "logical", test_function = is.logical)
  # Check numeric parameters.
  numeric_list <- list("pval_cutoff" = pval_cutoff,
                       "FC_cutoff" = FC_cutoff,
                       "pt.size" = pt.size,
                       "border.size" = border.size,
                       "font.size" = font.size,
                       "line_size" = line_size,
                       "n_genes" = n_genes)
  check_type(parameters = numeric_list, required_type = "numeric", test_function = is.numeric)
  # Check character parameters.
  character_list <- list("border.color" = border.color,
                         "font.type" = font.type,
                         "line_color" = line_color,
                         "plot.title" = plot.title,
                         "plot.subtitle" = plot.subtitle,
                         "plot.caption" = plot.caption,
                         "add_tags_for" = add_tags_for,
                         "order_tags_by" = order_tags_by,
                         "colors.use" = colors.use,
                         "plot.title.face" = plot.title.face,
                         "plot.subtitle.face" = plot.subtitle.face,
                         "plot.caption.face" = plot.caption.face,
                         "axis.title.face" = axis.title.face,
                         "axis.text.face" = axis.text.face,
                         "legend.title.face" = legend.title.face,
                         "legend.text.face" = legend.text.face)
  check_type(parameters = character_list, required_type = "character", test_function = is.character)

  check_colors(border.color, parameter_name = "border.color")
  check_colors(line_color, parameter_name = "line_color")
  check_colors(colors.use, parameter_name = "colors.use")

  check_parameters(parameter = font.type, parameter_name = "font.type")
  check_parameters(plot.title.face, parameter_name = "plot.title.face")
  check_parameters(plot.subtitle.face, parameter_name = "plot.subtitle.face")
  check_parameters(plot.caption.face, parameter_name = "plot.caption.face")
  check_parameters(axis.title.face, parameter_name = "axis.title.face")
  check_parameters(axis.text.face, parameter_name = "axis.text.face")
  check_parameters(legend.title.face, parameter_name = "legend.title.face")
  check_parameters(legend.text.face, parameter_name = "legend.text.face")


  assertthat::assert_that(order_tags_by %in% c("both", "pvalue", "logfc"),
                          msg = "Please use either both, pvalue or logfc in order_tags_by.")

  assertthat::assert_that(add_tags_for %in% c("both", "positive", "negative"),
                          msg = "Please use either both, positive or negative in add_tags_for")

  `%>%` <- magrittr::`%>%`
  colors <- do_ColorPalette(colors.use, tetradic = TRUE)
  names(colors) <- c("A", "C", "B", "D")

  if (!("gene" %in% colnames(de_genes))){
    data <- de_genes %>%
      tibble::rownames_to_column(var = "gene")
  } else {
    data <- de_genes
  }


  data <- data %>%
    tibble::as_tibble() %>%
    dplyr::select(c("p_val_adj", "avg_log2FC", "gene")) %>%
    dplyr::mutate("p_val_adj" = replace(.data$p_val_adj, .data$p_val_adj == 0, .Machine$double.xmin)) %>%
    dplyr::mutate(log_p = -log10(.data$p_val_adj)) %>%
    dplyr::select(-"p_val_adj")

  pval_cutoff <- -log10(pval_cutoff)
  data$color <- NA
  data$color[abs(data$avg_log2FC) >= FC_cutoff & data$log_p >= pval_cutoff] <- "A"
  data$color[abs(data$avg_log2FC) < FC_cutoff & data$log_p >= pval_cutoff] <- "B"
  data$color[abs(data$avg_log2FC) < FC_cutoff & data$log_p < pval_cutoff] <- "C"
  data$color[abs(data$avg_log2FC) >= FC_cutoff & data$log_p < pval_cutoff] <- "D"

  max_value <- max(abs(c(min(data$avg_log2FC), max(data$avg_log2FC))))
  x_lims <- c(-max_value, max_value)

  # Shuffle the data.
  data <- data[sample(rownames(data), nrow(data)), ]
  p <- data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$avg_log2FC,
                                           y = .data$log_p)) +
    ggplot2::geom_point(size = pt.size * border.size,
                        color = border.color) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = .data$color),
                        size = pt.size) +
    ggplot2::labs(title = plot.title,
                  subtitle = plot.subtitle,
                  caption = plot.caption) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4),
                                                  title.position = "top",
                                                  title.hjust = 0.5)) +
    ggplot2::xlim(x_lims) +
    ggplot2::xlab(expression(bold(paste("Avg. ", log["2"], "(FC)")))) +
    ggplot2::ylab(expression(bold(paste("-", log["10"], "(p-value adjusted)"))))

  if (isTRUE(plot_lines)){
    p <- p +
      ggplot2::geom_hline(yintercept = pval_cutoff,
                          color = line_color,
                          linewidth = line_size,
                          linetype = "dashed") +
      ggplot2::geom_vline(xintercept = FC_cutoff,
                          color = line_color,
                          linewidth = line_size,
                          linetype = "dashed") +
      ggplot2::geom_vline(xintercept = -FC_cutoff,
                          color = line_color,
                          linewidth = line_size,
                          linetype = "dashed")
  }

  if (isTRUE(add_gene_tags)){
    if (order_tags_by == "both"){
      data.up <- data %>%
        dplyr::arrange(dplyr::desc(.data$log_p),
                       dplyr::desc(.data$avg_log2FC)) %>%
        as.data.frame() %>%
        utils::head(n_genes)

      data.down <- data %>%
        dplyr::arrange(dplyr::desc(.data$log_p),
                       .data$avg_log2FC) %>%
        as.data.frame() %>%
        utils::head(n_genes)
    } else if (order_tags_by == "pvalue"){
      data.up <- data %>%
        dplyr::filter(.data$avg_log2FC > 0) %>%
        dplyr::arrange(dplyr::desc(.data$log_p),
                       dplyr::desc(.data$avg_log2FC)) %>%
        as.data.frame() %>%
        utils::head(n_genes)

      data.down <- data %>%
        dplyr::filter(.data$avg_log2FC < 0) %>%
        dplyr::arrange(dplyr::desc(.data$log_p)) %>%
        as.data.frame() %>%
        utils::head(n_genes)
    } else if (order_tags_by == "logfc"){
      data.up <- data %>%
        dplyr::arrange(dplyr::desc(.data$avg_log2FC)) %>%
        as.data.frame() %>%
        utils::head(n_genes)

      data.down <- data %>%
        dplyr::arrange(.data$avg_log2FC) %>%
        as.data.frame() %>%
        utils::head(n_genes)
    }

    if (add_tags_for == "both") {

      data.label <- dplyr::bind_rows(data.up, data.down)

    } else if (add_tags_for == "positive") {

      data.label <- data.up

    } else if (add_tags_for == "negative") {

      data.label <- data.down

    }

    if (base::isFALSE(use_labels)){
      p <- p +
        ggrepel::geom_text_repel(data = data.label,
                                 mapping = ggplot2::aes(label = .data$gene),
                                 max.overlaps = 1000,
                                 color = "black",
                                 fontface = "bold")
    } else if (isTRUE(use_labels)){
      p <- p +
        ggrepel::geom_label_repel(data = data.label,
                                  mapping = ggplot2::aes(label = .data$gene),
                                  max.overlaps = 1000,
                                  color = "black",
                                  fontface = "bold")
    }

  }
  p <- p +
    ggplot2::theme_minimal(base_size = font.size) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
                   plot.title = ggplot2::element_text(face = plot.title.face, hjust = 0),
                   plot.subtitle = ggplot2::element_text(face = plot.subtitle.face, hjust = 0),
                   plot.caption = ggplot2::element_text(face = plot.caption.face, hjust = 1),
                   legend.text = ggplot2::element_text(face = legend.text.face),
                   legend.title = ggplot2::element_text(face = legend.title.face),
                   panel.grid = ggplot2::element_blank(),
                   plot.title.position = "plot",
                   plot.caption.position = "plot",
                   text = ggplot2::element_text(family = font.type),
                   legend.position = "none",
                   legend.justification = "center",
                   axis.title.x = ggplot2::element_text(face = axis.title.face, color = "black"),
                   axis.title.y = ggplot2::element_text(face = axis.title.face, angle = 90, color = "black"),
                   axis.text = ggplot2::element_text(face = axis.text.face, color = "black"),
                   axis.line = ggplot2::element_line(color = "black"),
                   axis.ticks = ggplot2::element_line(color = "black"),
                   plot.background = ggplot2::element_rect(fill = "white", color = "white"),
                   panel.background = ggplot2::element_rect(fill = "white", color = "white"),
                   legend.background = ggplot2::element_rect(fill = "white", color = "white"))
  return(p)
}


## --- allPlot -------------------------------------------------------------- ##

allPlot <- function(object,
                    markers,
                    meta,
                    cluster,
                    split,
                    filename = "",
                    h = 7,
                    w = 15,
                    return = F,
                    save = T) {
  # --- Documentation ------------------------------------------------------- ##

  #' Documentation
  #' @title allPlot
  #'
  #' @description allPlot generates a dimplot, bar plot and volcano plot
  #'
  #' @param object Seurat object
  #' @param markers Dataframe. Output of FindMarkers in Seurat
  #' @param meta Character string. Metadata to use
  #' @param cluster Character string. Name of cluster to highlight
  #' @param split Character string. Name of split for barplot
  #' @param filename Characer string. Name of save directory location
  #' @param h Numeric. Height of saved file. Default is 7
  #' @param w Numeric. Width of saved file. Default is 15
  #' @param return Logical. Whether to return the plot or not. Default is False
  #' @param save Logical. Whether to save the plot or not. Default is True.
  #'
  #' @returns A ggplot grid
  #'
  #' @import SCpubr
  #' @import ggplot2
  #' @importFrom cowplot plotgrid
  #'
  #' @export

  # --- End Documentation --------------------------------------------------- ##


  # Get metadata from object and convert to data frame and order by frequency
  data <- as.data.frame(table(object@meta.data[meta]))

  data <- data[order(data$Freq, decreasing = T),]

  # Set color scheme and set organ of choice color as steelblue
  color <- do_ColorPalette("steelblue", nrow(data))

  rest <- setdiff(data[,1], cluster)

  names(color) <- c(cluster, rest)

  my_levels <- data[,1]

  # Change factor levels of active ident
  object <- SetIdent(object, value = meta)

  object@meta.data[[meta]] <- factor(x = object@meta.data[[meta]], levels = my_levels)

  leg_n <- nrow(data)

  leg_row <- ceiling(leg_n/3)

  # DimPlot of object
  dPlot <- do_DimPlot(sample = object,
                      group.by = meta,
                      legend.icon.size = 5,
                      font.size = 20,
                      colors.use = color,
                      legend.ncol = 3,
                      legend.nrow = leg_row)

  # Percentage stacked bar chart showing population percentages
  bPlot <- do_BarPlot(object,
                      group.by = meta,
                      split.by = split,
                      position = "fill",
                      flip = FALSE,
                      legend.position = "none",
                      colors.use = color,
                      font.size = 20) + theme(axis.text.x = element_blank(),
                                              axis.title.x = element_blank())

  # Volcano plot of DE genes
  vPlot <- do_VolcanoPlot_modified(sample = object,
                          de_genes = markers,
                          FC_cutoff = 0.25,
                          colors.use = c("darkred"),
                          n_genes = 5,
                          add_gene_tags = T,
                          add_tags_for = 'positive',
                          use_labels = T,
                          font.size = 20,
                          pval_cutoff = 1)

  # Use plot_grid to combine all plots
  p <- plot_grid(dPlot, bPlot, vPlot, nrow=1, ncol=3, rel_widths = c(4,1,4))

  if (save == T) {

    ggsave(plot = p,
           filename = filename,
           width = w,
           height = h)

  }


  if (return == T) {

    return(p)

  }

}


# Citation Christian Blötner, diffcor,
# https://github.com/cran/diffcor/blob/master/R/diffcor.dep.R

diffcor.dep <- function(r12, r13, r23, n, cor.names = NULL,
                        alternative = c("one.sided", "two.sided"), digit = 3) {

  # --- Documentation ------------------------------------------------------- ##

  #' @title diffcor.dep
  #'
  #' @description Test of differences between two dependent correlations
  #'
  #' @param r12 empirically observed correlation between first and second construct
  #' @param r13 empirically observed correlation between first and third construct
  #' @param r23 empirically observed correlation between second and third construct
  #' @param n sample size the correlations are based on
  #'
  #' @return Fisher's z-value and corresponding p-values
  #'
  #' @export

  z12 <- atanh(r12)
  z13 <- atanh(r13)
  r1 <- (r12 + r13) / 2

  Cov.dep.a <- 1 / ((1 - (r1 ^ 2))^2)
  Cov.dep.b <- r23 * (1 - (2 * (r1^2)))
  Cov.dep.c <- .50 * (r1^2)
  Cov.dep.d <- 1 - (2 * (r1^2)) - (r23^2)
  Cov.dep <- (Cov.dep.a * Cov.dep.b) - (Cov.dep.c * Cov.dep.d)
  SE.dep <- sqrt((2 - (2 * Cov.dep)) / (n - 3))
  diff.z.dep <- round(((z12 - z13) / SE.dep), digit)

  alternative <- match.arg(alternative);
  tside <- switch(alternative, one.sided = 1, two.sided = 2)

  if (tside == 1){
    p <- format(round(1 - pnorm(abs(diff.z.dep)), digit), scientific = F)}

  if (tside == 2){
    p <- format(round(2 * pnorm(-abs(diff.z.dep)), digit), scientific = F)}

  res.dep <- data.frame(round(r12, digit), round(r13, digit), round(r23, digit),
                        diff.z.dep, p)
  rownames(res.dep) <- cor.names
  colnames(res.dep) <- c("r12", "r13", "r23", "z", "p")

  return(res.dep)
}

## --- fisherZtransform ----------------------------------------------------- ##

# Function to perform the test
# Website https://www.ibm.com/support/pages/differences-between-correlations
# Citations :
# Sheshkin, D.J. (2004). Handbook of Parametric and Nonparametric Statistical Procedures (3rd Ed.). Boca Raton FL: Chapman & Hall/CRC.
# Zar, J.H. (1999). Biostatistical Analysis(4th Ed.). Upper Saddle River NJ: Prentice Hall.

fisherZtransform <- function(r1, r2, n1, n2, method = 'spearman') {

  # --- Documentation ------------------------------------------------------- ##

  #' @title fisherZtransform
  #'
  #' @description This is a function to conduct a fisher Z transformation on Pearson's or Spearman's Rho correlation
  #'
  #' @param r1,r2 Float. Correlations to transform
  #' @param n1,n2 Numeric. Number of samples
  #' @param method Character string. Name of correlation statitic used, Spearman, Pearson
  #'
  #' @export

  # --- End of Documentation ------------------------------------------------ ##


  if (method == 'spearman') {

    cor_factor <- 1.06 # Correction factor for spearman Sheshkin, 2004; Zar, 1999

  } else if (method == 'pearson') {

    cor_factor <- 1
  }
  # Function to calculate Fisher Z transformation
  fisher_z <- function(rho) {
    0.5 * log((1 + rho) / (1 - rho))
  }

  # Function to calculate standard error of Fisher Z difference
  se_diff <- function(r1, r2, n1, n2) {
    sqrt((cor_factor / (n1 - 3)) + (cor_factor / (n2 - 3)))
  }

  z1 <- fisher_z(r1)

  z2 <- fisher_z(r2)

  zdiff <- z1 - z2

  se_diff_val <- se_diff(r1, r2, n1, n2)

  ztest <- (z1 - z2) / se_diff_val

  alpha <- 2 * (1 - pnorm(abs(ztest), 0, 1))

  return(list(diff = zdiff, test = ztest, p = alpha))
}


## --- doSpatialFeaturPlot -------------------------------------------------- ##

diffCorrPlot_n <- function(list,
                           features,
                           ident.1 = 'endocrine',
                           ident.2 = 'exocrine',
                           assay = 'SCT',
                           slot = 'data',
                           stat = 'spearman',
                           plot = 'bar',
                           custom_group = F,
                           group_var = "",
                           group_pal = "",
                           order = F,
                           flip = T) {

  # --- Documentation ------------------------------------------------------- ##

  #' @title diffCorrPlot_n
  #'
  #' @description This is a function that outputs a bar graph or box plot to demonstrate differences in correlation across multiple seperate objects for two unrelated identities
  #'
  #' @param list List of Seurat Objects
  #' @param fatures List of character strings. Features to plot
  #' @param ident.1,ident.2  Character strings. Feature to conduct correlation on
  #' @param assay Character string. Which assay to pull data from if data are expression values
  #' @param slot Character string. Which slot to pull data from if data are expression values.
  #' @param stat Character string. Name of correlation statistic used, Spearman, Pearson
  #' @param plot Character string. Type of plot to output, "bar" for bar graph of Z scores or "box" for box plots of correlations. Default is "bar'.
  #' @param custom_group Logical. Whether to include custom grouping for plots, only for bar graph. Default is F
  #' @param group_var Dataframe. feature-group dataframe for custom grouping.
  #' @param flip Logical. If T. Then the coordinates are flipped
  #' @export
  #'
  #' @import Seurat
  #' @import dplyr
  #' @import tidyr
  #'

  # --- End of Documentation ------------------------------------------------ ##


  object_list <- list

  # Set up empty dataframes to add to

  corr_df_all <- data.frame(row.names = features, "feature" = features)

  corr_df_z_all <- data.frame(row.names = features, "feature" = features)

  for (i in 1:length(object_list)) {

    obj <- object_list[[i]]

    meta <- FetchData(obj, vars = c(features, ident.1, ident.2))

    cor1 <- cor(meta[,ident.1], meta[,features], method = stat)

    cor2 <- cor(meta[,ident.2], meta[,features], method = stat)

    corr_df <- rbind(cor1, cor2) %>%
      t() %>%
      as.data.frame()

    colnames(corr_df) <- c(paste0("N",i,"_",ident.1),
                           paste0("N",i,"_",ident.2))

    corr_df$feature <- features

    corr_df_all <- corr_df_all %>%
      left_join(corr_df, by = "feature")

    r1 <- corr_df[,1]

    r2 <- corr_df[,2]

    n1 <- n2 <- nrow(meta)

    corr_df_z<- data.frame(row.names = features, "feature" = features)

    corr_df_z[,paste0("N",i,"_Z")] <- fisherZtransform(r1, r2, n1, n2, method = stat)$test

    corr_df_z_all <- corr_df_z_all %>%
      left_join(corr_df_z, by = "feature")

  }

  # Pivot longer
  corr_df_longer <- pivot_longer(corr_df_all,
                                 cols = starts_with('N')) %>%
    mutate(feature = factor(feature, levels = corr_df_all$feature))

  corr_df_longer[,c("N", "ident")] <- str_split_fixed(corr_df_longer$name,
                                                      pattern = "_",
                                                      n = 2)

  p_box <- ggboxplot(corr_df_longer,
                     x = "feature",
                     y = "value",
                     fill = "ident",
                     line.color = "black",
                     line.size = 0.4,
                     palette = c("#1F78B4FF", "#33A02CFF")) +
    theme( plot.title = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_text(size=8, face = "bold", colour = "black"),
           axis.text.x = element_text(size =8, angle = 90, face = "bold", colour = "black", vjust = 0.5),
           axis.text.y = element_text(size =8, face = "bold", colour = "black"),
           legend.title = element_text(size=8, face = "bold", colour = "black", angle = 0),
           legend.text = element_text(size=8, face = "bold", colour = "black"),
           legend.title.align = 0.5,
           legend.position = 'none') +
    labs(y = paste(str_to_title(stat), 'Correlation', sep = " "))



  corr_df_z_longer <- pivot_longer(corr_df_z_all, cols = starts_with('N')) %>%
    mutate(feature = factor(feature, levels = corr_df_z_all$feature))

  corr_df_z_longer[,c("N")] <- str_split_fixed(corr_df_z_longer$name, pattern = "_", n = 2)[,1]

  if (length(features) < 6) {

    pal_bar <- c("#E31A1CFF", "#FF7F00FF", "#33A02CFF","#1F78B4FF", "#6A3D9AFF")

  } else {

    pal_orig<- c("#E31A1CFF", "#FF7F00FF", "#33A02CFF","#1F78B4FF", "#6A3D9AFF")

    pal_bar <- c(pal_orig,
                 setdiff(pal_orig,
                         do_ColorPalette("#E31A1CFF", n = length(features))))
  }

  if (custom_group == T) {

    corr_df_z_longer <- corr_df_z_longer %>%
      left_join(group_var, by = "feature")

    if (order == T) {

      x_order <- corr_df_z_longer %>%
        filter(!is.na(value)) %>%
        group_by(feature) %>%
        summarize(mean = mean(value)) %>%
        arrange(mean) %>%
        pull(feature)

    } else {
      x_order <- NULL
    }

    if (flip == T) {
      p_bar <- ggbarplot(corr_df_z_longer,
                         x = "feature",
                         y = "value",
                         fill = "group",
                         add = "mean_se",
                         palette = group_pal,
                         order = x_order) +
        theme( plot.title = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_text(size=8, face = "bold", colour = "black"),
               axis.text.x = element_text(size =8, angle = 45, face = "bold", colour = "black", vjust = 0.5),
               axis.text.y = element_text(size =8, face = "bold", colour = "black"),
               legend.title = element_text(size=8, face = "bold", colour = "black", angle = 0),
               legend.text = element_text(size=8, face = "bold", colour = "black"),
               legend.title.align = 0.5,
               legend.position = 'right') + coord_flip() +
      labs(y = "Fisher Z score")

    } else if (flip == F) {

      p_bar <- ggbarplot(corr_df_z_longer,
                         x = "feature",
                         y = "value",
                         fill = "group",
                         add = "mean_se",
                         palette = group_pal,
                         order = x_order) +
        theme( plot.title = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_text(size=8, face = "bold", colour = "black"),
               axis.text.x = element_text(size =8, angle = 45, face = "bold", colour = "black", vjust = 0.5),
               axis.text.y = element_text(size =8, face = "bold", colour = "black"),
               legend.title = element_text(size=8, face = "bold", colour = "black", angle = 0),
               legend.text = element_text(size=8, face = "bold", colour = "black"),
               legend.title.align = 0.5,
               legend.position = 'right') +
      labs(y = "Fisher Z score")

    }

  } else {

    p_bar <- ggboxplot(corr_df_z_longer,
                       x = "feature",
                       y = "value",
                       fill = "feature",
                       color = 'black',
                       palette = pal_bar) +
      scale_fill_manual(values = pal_bar) +
      theme( plot.title = element_blank(),
             axis.title.x = element_blank(),
             axis.title.y = element_text(size=8, face = "bold", colour = "black"),
             axis.text.x = element_text(size =8, angle = 45, face = "bold", colour = "black", vjust = 0.5),
             axis.text.y = element_text(size =8, face = "bold", colour = "black"),
             legend.title = element_text(size=8, face = "bold", colour = "black", angle = 0),
             legend.text = element_text(size=8, face = "bold", colour = "black"),
             legend.title.align = 0.5,
             legend.position = 'none') +
      labs(y = "Fisher Z score")

  }



  if (plot == 'bar') {

    return(list(corr_df_z_longer, p_bar))

    print()

  } else if (plot == 'box') {

    return(list(corr_df_longer, p_box))
  }

}



## --- doSpatialFeaturPlot -------------------------------------------------- ##

min_euc_dist <- function(object,
                         metadata = 'region',
                         meta_ident = 'endocrine',
                         add_meta_name = paste0('euc_dist_',
                                                meta, "_",
                                                meta_ident)) {

  # --- Documentation ------------------------------------------------------- ##

  #' @title min_euc_dist
  #'
  #' @description This is a function to annotate the minimum distance from specified voxels
  #'
  #' @param object Seurat spatial object.
  #' @param meta  Character string. Which meta data to use for specifying voxels
  #' @param meta_ident Character string. Which meta_identity to use for specifying voxels
  #'
  #' @export

  object@meta.data$row <- object@images$slice1@coordinates$row

  object@meta.data$col <- object@images$slice1@coordinates$col

  meta <- FetchData(object, vars = c(metadata, "row", "col"))

  colnames(meta)[1] <- 'feature'

  meta_point <- meta %>%
    filter(feature == meta_ident) %>%
    select(-feature)

  # Create an empty vector to store the distances
  nearest_dist <- c()

  # Loop through each sample in the dataframe
  for (i in 1:nrow(meta)) {

    sample <- meta[i, c("row", "col")]

    point_coord <- meta_point

    point_coord[,"row"] <- meta_point[,"row"] - sample[,"row"]

    point_coord[,"col"] <- meta_point[,"col"] - sample[,"col"]

    # Calculate distances to all "Endocrine" samples
    dist <- sqrt(rowSums((point_coord)^2))

    # Find the minimum distance
    min_dist <- min(dist)

    # Append the minimum distance to the vector
    nearest_dist <- c(nearest_dist, min_dist)
  }

  object@meta.data[[add_meta_name]] <- nearest_dist

  return(object)

}
