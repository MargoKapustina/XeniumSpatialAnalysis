#' Xenium 1-D UMAP Embeddings Plotted in Situ
#' @author Margarita Kapustina
#'
#' @description Plots: 1-dimensional UMAP, the corresponding histogram of UMAP_1 embedding values, and
#'    the UMAP_1 embedding values in situ within specified FOV.
#'
#' @param object a Xenium object
#' @param FOV FOV to plot UMAP_1 embedding values in
#' @param EmbeddingsPlotTitle title for UMAP_1 Embeddings plot
#' @param HistogramPlotTitle title for UMAP_1 Embeddings Histogram plot
#' @param inSituPlotTitle title for UMAP_1 Embeddings in situ plot
#' @param save_plot Option to save plot as .eps in working directory (TRUE, FALSE)
#'
#' @return Dataframe containing Cell IDs, coordinates (X,Y) and computed UMAP_1 embeddings values for cells in specified FOV.
#' #Plots 1D UMAP, corresponding histogram of UMAP_1 embeddings values, and UMAP_1 embeddings values plotted _in situ_.
#'
#' @import ggplot2
#' @import gridExtra
#' @import viridis
#' @import Seurat
#' @import SeuratObject
#' @import tibble
#'
#' @examples
#' AD_df = plotUMAP1_inSitu(object = AD, FOV = 'fov')
#' @export

plotUMAP1_inSitu <- function(object,
                             FOV = 'X1fov',
                             EmbeddingsPlotTitle = 'UMAP_1 embedding for all\n cells in subset',
                             HistogramPlotTitle = 'Distribution of UMAP_1 embedding\n values all cells in subset',
                             inSituPlotTitle = 'UMAP_1 embedding in situ\n for representative FOV',
                             save_plot = FALSE){

  #get current default assay to then set back later
  current_default_assay = SeuratObject::DefaultAssay(object)

  # Get UMAP_1 embeddings from Xenium object
  cell.embeddings <- object@reductions$umap@cell.embeddings
  cell.embeddings.df <- as.data.frame(cell.embeddings)
  names(cell.embeddings.df)[1] <- "UMAP_1"

  # Add Cell_ID as a column
  cell.embeddings.df <- tibble::rownames_to_column(cell.embeddings.df, "Cell_ID")
  #rename cell.embeddings.df as df for plotting use
  df = cell.embeddings.df
  # Add a column 'y' with default value 0
  df$y <- 0

  # Plot continuum of UMAP_1 embedding values
  embeddings_plot <- ggplot2::ggplot(df, aes(x = UMAP_1, y = y, colour = UMAP_1)) +
    geom_point(size = 0.5) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_color_viridis_c() +
    geom_jitter() +
    ggtitle(EmbeddingsPlotTitle)

  #histogram of embeddings across all slices
  embeddings_histogram = ggplot2::ggplot(df, aes(x = UMAP_1, fill = ..x..)) +
    geom_histogram(bins = 40, colour = 'blue')+
    ggtitle(HistogramPlotTitle)+
    viridis::scale_fill_viridis()
  #set default FOV to fetch coordinates for FOV of choice
  SeuratObject::DefaultFOV(object) <- FOV

  # Get tissue coordinates
  df <- Seurat::GetTissueCoordinates(object)

  # Merge tissue coordinates (FOV of interest only) with the dataframe containing 1D UMAP embeddings
  cells_to_get <- df[["cell"]]
  cell.embeddings.df_subset <- cell.embeddings.df[cell.embeddings.df$Cell_ID %in% cells_to_get, ]
  names(cell.embeddings.df_subset)[1] <- 'cell'
  coordsEmbeddingsdf <- merge(df, cell.embeddings.df_subset, by = 'cell')

  # Plot embedding values by colour in FOV of choice
  inSituPlot = ggplot2::ggplot(coordsEmbeddingsdf, aes(x = x, y = y)) +
    geom_point(aes(colour = UMAP_1), size = 1.5) +
    scale_color_viridis_c()+
    ggtitle(inSituPlotTitle) +
    coord_fixed()
  # Print the plots side by side
  # print(embeddings_plot|embeddings_histogram|inSituPlot)
  print("Generating plots...")
  gridExtra::grid.arrange(
    gridExtra::arrangeGrob(embeddings_plot, embeddings_histogram, ncol = 1),
    inSituPlot,
    ncol = 2,
    widths = c(0.5, 0.5)
  )
  # saving plot
  if (save_plot) {
    print('Saving plot...')
    dev.copy2eps(file = 'UMAP1EmbeddingsInSitu_plot.eps')
    print("Plot saved to working directory.")}

  #change default assay back to whatever it was before running function
  SeuratObject::DefaultAssay(object) <- current_default_assay

  # Return the merged data frame
  return(coordsEmbeddingsdf)
}

