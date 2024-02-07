#' Xenium 1-D UMAP Embeddings Plotted in Situ
#' @author Margarita Kapustina
#'
#' @description Plots: 1-dimensional UMAP, the corresponding histogram of UMAP_1 embedding values,
#'    the UMAP_1 embedding values in situ within specified FOV, and computes the midline intersecting cells within specified FOV.
#'    note: please compute UMAP_1 embeddings beforehand.
#'
#' @param object a Xenium object
#' @param FOV FOV to plot UMAP_1 embedding values in
#' @param degs User-defined angle, defined in degrees, that intersects centrepoint in midline computation
#' @param EmbeddingsPlotTitle title for UMAP_1 Embeddings plot
#' @param HistogramPlotTitle title for UMAP_1 Embeddings Histogram plot
#' @param inSituPlotTitle title for UMAP_1 Embeddings in situ plot
#' @param save_plot Option to save plot as .eps in working directory (TRUE, FALSE)
#'
#' @return Dataframe containing Cell IDs, coordinates (X,Y) and computed UMAP_1 embeddings values for cells in specified FOV.
#' #Plots 1D UMAP, corresponding histogram of UMAP_1 embeddings values, and UMAP_1 embeddings values plotted in situ.
#'
#' @import ggplot2
#' @import gridExtra
#' @import viridis
#' @import Seurat
#' @import SeuratObject
#' @import tibble
#'
#' @examples
#' AD_df = plotUMAP1_inSituMidline(object = AD, FOV = 'fov', degs = 65)
#' @export

plotUMAP1_inSituMidline <- function(object,
                                    FOV = 'X1fov',
                                    degs = 45,
                                    save_plot = FALSE,
                                    EmbeddingsPlotTitle = 'UMAP_1 embedding for all\n cells in subset',
                                    HistogramPlotTitle = 'Distribution of UMAP_1 embedding\n values all cells in subset',
                                    inSituPlotTitle = 'UMAP_1 embedding in situ\n for representative FOV'){
  #store .getCentre() Function
  .getCentre <- function(df,doMean=F){
    if(doMean){
      dX <- mean(df$X)
      dY <- mean(df$Y)
    }else{
      dX <- median(df$X)
      dY <- median(df$Y)
    }

    dX <- df$X-dX
    dY <- df$Y-dY

    df <- cbind(df,dX,dY)

    return(df)
  }

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
  embeddings_plot <- ggplot2::ggplot(df, ggplot2::aes(x = UMAP_1, y = y, colour = UMAP_1)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::geom_jitter() +
    ggplot2::ggtitle(EmbeddingsPlotTitle)

  #histogram of embeddings across all slices
  embeddings_histogram = ggplot2::ggplot(df, ggplot2::aes(x = UMAP_1, fill = ggplot2::after_stat(x))) +
    ggplot2::geom_histogram(bins = 40, colour = 'blue')+
    ggplot2::ggtitle(HistogramPlotTitle)+
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
  inSituPlot = ggplot2::ggplot(coordsEmbeddingsdf, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(colour = UMAP_1), size = 1.5) +
    ggplot2::scale_color_viridis_c()+
    ggplot2::ggtitle(inSituPlotTitle) +
    ggplot2::coord_fixed()

  #plot UMAP_1 embedding with nucleus midline as well######################
  names(coordsEmbeddingsdf)[names(coordsEmbeddingsdf) == 'x'] <- 'X'
  names(coordsEmbeddingsdf)[names(coordsEmbeddingsdf) == 'y'] <- 'Y'

  #compute the centre of the coordinates
  centre = coordsEmbeddingsdf %>% .getCentre(doMean = TRUE)
  #set slope of line
  degs = degs
  rads <- degs*pi/180
  m <- tan(rads)
  #run distance operation and make new df "myplot"
  dLine <- (-m*centre$dX+centre$dY)/sqrt(1+m^2)
  myplot <- cbind(centre,dLine)
  #plot out UMAP_1 embedding with nucleus midline
  gg <- ggplot2::ggplot(myplot,ggplot2::aes(x=dX,y=dY))
  gg <- gg + ggplot2::geom_point(ggplot2::aes(colour= UMAP_1), size = 1.5)
  gg <- gg + ggplot2::geom_abline(slope=m)
  gg <- gg + ggplot2::coord_fixed() +ggplot2::ggtitle(inSituPlotTitle)+ggplot2::scale_color_viridis_c()
  print(gg)

  # Print the plots side by side
  # print(embeddings_plot|embeddings_histogram|gg)
  print("Generating plots...")
  print("If in situ spatial midline looks off, please adjust 'degs' parameter!")
  gridExtra::grid.arrange(
    gridExtra::arrangeGrob(embeddings_plot, embeddings_histogram, ncol = 1),
    gg,
    ncol = 2,
    widths = c(0.5, 0.5)
  )
  # saving plot
  if (save_plot) {
    print('Saving plot...')
    dev.copy2eps(file = 'UMAP1EmbeddingsInSituMidline_plot.eps')
    print("Plot saved to working directory.")}

  #change default assay back to whatever it was before running function
  SeuratObject::DefaultAssay(object) <- current_default_assay

  # Return the merged data frame
  return(coordsEmbeddingsdf)
}
