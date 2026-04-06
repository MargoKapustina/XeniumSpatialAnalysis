#' Plot Cluster Frequencies and Proportions Stratified by Metadata
#'
#' Creates two bar plots showing (1) raw cell counts and (2) normalized proportions
#' of clusters across a metadata variable. Optionally facets by a third metadata variable.
#'
#' @param obj A Seurat object with metadata variables
#' @param cluster_var Character. Metadata variable containing cluster labels
#' @param group_var Character. Metadata variable used for grouping (fill variable)
#' @param facet_var Optional character. Metadata variable to facet plots by (Default is NULL)
#' @param colors Character vector of colors. Used if number of groups is <= length(colors)
#'
#' @return A list with two ggplot objects:
#' \item{raw}{Bar plot showing raw cell counts.}
#' \item{norm}{Bar plot showing proportions.}
#'
#' @import ggplot2
#' @import Seurat
#' 
#' @examples
#' #generate plots
#' plots <- plotMetadataFrequency(obj = myObj, cluster_var = "seurat_clusters", group_var = "condition", facet_var = "sex")
#' #print plots
#' plots$raw + plots$norm
#'
#' @export

plotMetadataFrequency <- function(obj, cluster_var, group_var, facet_var = NULL, colors = c("magenta", "grey")) {
  
  # build dataframe
  df <- data.frame(
    cluster = obj[[cluster_var]][,1],
    group   = obj[[group_var]][,1]
  )
  
  # optional facet
  if (!is.null(facet_var)) {
    df$facet <- obj[[facet_var]][,1]
  }
  
  # factors
  df$cluster <- factor(df$cluster, levels = sort(unique(df$cluster)))
  df$group <- factor(df$group)
  
  n_groups <- length(levels(df$group))
  
  # choose palette
  if (n_groups <= length(colors)) {
    fill_scale <- ggplot2::scale_fill_manual(values = colors[1:n_groups])
  } else {
    fill_scale <- ggplot2::scale_fill_brewer(palette = "Set2")
  }
  
  # RAW plot
  p_raw <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, fill = group)) +
    ggplot2::geom_bar(position = "dodge") +
    fill_scale +
    ggplot2::labs(
      x = cluster_var,
      y = "Cell Count",
      fill = group_var,
      title = paste("Cell Counts by", cluster_var, "and\n", group_var)
    ) +
    ggplot2::theme_minimal()
  
  # NORMALIZED plot
  p_norm <- ggplot2::ggplot(df, ggplot2::aes(x = cluster, fill = group)) +
    ggplot2::geom_bar(position = "fill") +
    fill_scale +
    ggplot2::labs(
      x = cluster_var,
      y = "Proportion",
      fill = group_var,
      title = paste("Proportion of", group_var, "within each\n", cluster_var)
    ) +
    ggplot2::theme_minimal()
  
  # facet if requested
  if (!is.null(facet_var)) {
    p_raw <- p_raw + ggplot2::facet_wrap(~facet)
    p_norm <- p_norm + ggplot2::facet_wrap(~facet)
  }
  
  return(list(raw = p_raw, norm = p_norm))
}
