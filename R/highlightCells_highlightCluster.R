#' Highlight Cluster
#' @author Margarita Kapustina
#'
#' @description Highlight cluster(s) of choice in a Xenium object, within specified FOV(s). Cluster of interest labelled "this cluster"
#'      and all other clusters labelled "all else". Supported for clusters within the same object.
#'
#' @param object a Xenium object
#' @param cluster_id Identity of cluster to highlight
#' @param size Size of the highlighted cells in the plot
#' @param alpha_value Alpha (transparency) value of the highlighted cells
#' @param FOV Character vector specifying the fields of view (FOVs) to use for plotting
#' @param color_palette Color palette for the plot
#' @param save_plot Option to save plot as .pdf in working directory (TRUE, FALSE)
#'
#'
#' @return An ImageDimPlot with highlighted cells.
#'
#' @examples
#' highlightCluster(obj, cluster_id = 1, FOV = c('X1fov', 'X2fov'))
#'
#' @import Seurat
#' @export
highlightCluster <- function(object, cluster_id, size = 3, alpha_value = 0.5, FOV = c("fov", 'X1fov'), color_palette = 'glasbey', save_plot = FALSE) {
  toPlot <- Seurat::WhichCells(object, idents = cluster_id)
  object@meta.data$ClusterToPlot <- ifelse(colnames(object) %in% toPlot, "this_cluster", "all_else")
  p1 = Seurat::ImageDimPlot(object, group.by = 'ClusterToPlot', size = size, alpha = alpha_value, fov = FOV, cols = color_palette)
  # saving plot
  if (save_plot) {
    print('Saving plot...')
    dev.copy2pdf(file = 'ClusterHighlighted.pdf')
    print("Plot saved to working directory.")}
  print(p1)
  #remove the ClusterToPlot column
  object@meta.data$ClusterToPlot = NULL
  return(p1)
}


#' Highlight Cells
#' @author Margarita Kapustina
#'
#' @description Highlight cells from a Xenium object of interest, within another Xenium object, within specified FOV(s).
#'    Cells from object of interest labelled as "these cells", other cells labelled as "all else". Suggested to use when subsetting Xenium objects.
#'
#' @param highlight_obj Xenium object that you want to highlight cells from, labelled: "these cells"
#' @param within_obj Xenium object that you want show highlighted cells within, labelled: "all else" (ie cells plotted but not highlighted)
#' @param size Size of the highlighted cells in the plot
#' @param alpha_value Alpha (transparency) value of the highlighted cells
#' @param FOV Character vector specifying the fields of view (FOVs) to use for plotting
#' @param color_palette Color palette for the plot
#' @param save_plot Option to save plot as .pdf in working directory (TRUE, FALSE)
#'
#' @return An ImageDimPlot with highlighted cells.
#'
#' @examples
#' highlightCells(highlight_obj = atn, within_obj = ExcitatoryNeurons, FOV = c('X1fov'))
#'
#' @import Seurat
#' @export

#highlight cells from one Xenium object, in another (ie good to use for subsets)
#object 1 is the cluster that you want to highlight - ie "these cells"
#object 2 will show you where its cells are plotted, but will not highlight them - ie "all else" cells.
highlightCells <- function(highlight_obj, within_obj, size = 3,
                           alpha_value = 0.5, FOV = c("fov", 'X1fov'),
                           color_palette = 'glasbey',
                           save_plot = FALSE) {
  toPlot <- Seurat::WhichCells(highlight_obj)
  within_obj@meta.data$CellsToPlot <- ifelse(colnames(within_obj) %in% toPlot, "these_cells", "all_else")
  p1 = Seurat::ImageDimPlot(within_obj, group.by = 'CellsToPlot', size = size, alpha = alpha_value, fov = FOV, cols = color_palette)
  # saving plot
  if (save_plot) {
    print('Saving plot...')
    dev.copy2pdf(file = 'ObjectHighlighted.pdf')
    print("Plot saved to working directory.")}
  print(p1)
  #remove the ClusterToPlot column
  within_obj@meta.data$CellsToPlot = NULL
  return(p1)
}
