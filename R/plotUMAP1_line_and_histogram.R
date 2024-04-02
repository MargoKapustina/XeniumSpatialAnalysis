#' Fetch UMAP_1 embeddings values for cells and their distance away from Spatial Midline
#' @author Margarita Kapustina
#'
#' @description Fetches UMAP_1 embeddinga values for cells, and their distance away from the computed spatial midline. Supports one FOV at a time.
#'
#' @param object a Xenium object
#' @param FOV FOV to extract coordinates and UMAP_1 embeddings values from
#' @param degs User-defined angle, defined in degrees, that intersects centrepoint in midline computation
#' @param binNumber Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth
#' @param binwidth Width of bins for gene expression averaging (lines). Suggested binwidth = 40 for 40micron bins
#'
#' @return Dataframe containing Cell IDs, coordinates (X,Y), UMAP_1 embedding values,
#'    and computed distance away from spatial midline of each cell, in specified FOV.
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import tibble
#'
#' @examples
#' UMAP1_midline_data_fov <- getUMAP1_MidlineData(AD, FOV= c('fov'), degs = 65)
#' UMAP1_midline_data_X1fov <- getUMAP1_MidlineData(AD, FOV= c('X1fov'), degs = 45)
#' @export

getUMAP1_MidlineData <- function(object, FOV, degs = 45, binNumber = 15, binwidth = 40) {

  # Fetch UMAP_1 embeddings for Xenium obj
  cell.embeddings = object@reductions$umap@cell.embeddings
  cell.embeddings.df = as.data.frame(cell.embeddings)
  names(cell.embeddings.df)[1] <- "UMAP_1"
  cell.embeddings.df <- tibble::rownames_to_column(cell.embeddings.df, "Cell_ID")

  #get current default assay to then set back later
  current_default_assay = SeuratObject::DefaultAssay(object)

  #set default FOV to fetch coordinates for FOV of choice
  SeuratObject::DefaultFOV(object) <- FOV

  # Get tissue coordinates
  df <- Seurat::GetTissueCoordinates(object)

  # Merge tissue coordinates with the dataframe containing 1D UMAP embeddings
  cells_to_get <- df[["cell"]]
  cell.embeddings.df_subset <- cell.embeddings.df[cell.embeddings.df$Cell_ID %in% cells_to_get, ]
  names(cell.embeddings.df_subset)[1] = 'cell'
  coordsEmbeddingsdf = merge(df, cell.embeddings.df_subset, by = 'cell')

  # Prepare data for plotting UMAP_1 embedding with nucleus midline
  names(coordsEmbeddingsdf)[names(coordsEmbeddingsdf) == 'x'] <- 'X'
  names(coordsEmbeddingsdf)[names(coordsEmbeddingsdf) == 'y'] <- 'Y'

  # Store .getCentre() Function
  .getCentre <- function(df, doMean = FALSE) {
    if(doMean) {
      dX <- mean(df$X)
      dY <- mean(df$Y)
    } else {
      dX <- median(df$X)
      dY <- median(df$Y)
    }

    dX <- df$X - dX
    dY <- df$Y - dY

    df <- cbind(df, dX, dY)

    return(df)
  }

  # Compute the centre of the coordinates
  centre = coordsEmbeddingsdf %>% .getCentre(doMean = TRUE)
  # Set slope of the line
  rads <- degs * pi / 180
  m <- tan(rads)
  # Run distance operation and make a new df "myplot"
  dLine <- (-m * centre$dX + centre$dY) / sqrt(1 + m^2)
  myplot <- cbind(centre, dLine)
  #plot out UMAP_1 embedding with nucleus midline
  gg <- ggplot2::ggplot(myplot,ggplot2::aes(x=dX,y=dY))
  gg <- gg + ggplot2::geom_point(ggplot2::aes(colour= UMAP_1), size = 1.5)
  gg <- gg + ggplot2::geom_abline(slope=m)
  gg <- gg + ggplot2::coord_fixed() +ggplot2::ggtitle('UMAP_1 values across Spatial Midline')+ggplot2::scale_color_viridis_c()
  print(gg)

  # Print the plots side by side
  # print(embeddings_plot|embeddings_histogram|gg)
  print("Generating plots...")
  print("If in situ midline looks off, please adjust 'degs' parameter!")

  # ggplot2 Visualization
  #show UMAP_1 embeddings as a line graph, as a distance away from midline
  p1 = ggplot2::ggplot(myplot,ggplot2::aes(x = dLine)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count / 100)), bins = binNumber, boundary= 0, colour = 'green')+
    ggplot2::coord_cartesian(xlim = c(-500, 500), ylim =c(-6,8))+
    ggplot2::xlab('Distance away from Spatial Midline')+
    ggplot2::ylab('Cells per bin/100')+
    ggplot2::stat_summary_bin(ggplot2::aes(x = dLine, y= UMAP_1),fun = mean, binwidth = binwidth, geom = 'line', colour = 'blue') +
    #stat_bin(bins = 7, boundary = 0,geom = 'text', aes(label = ..count../100), vjust = 0, size = 3, color = 'black')+
    ggplot2::ggtitle('UMAP_1 embeddings\nbins are mean number of cells\nper bin /100 in this FOV ONLY')
  print(p1|gg)

  # Calculate the number of cells in each bin
  bin_counts <- ggplot2::ggplot_build(p1)$data[[1]]
  cat("Number of cells in each bin (note: plotted are counts/100):
Beginning on left-most bin...\n")
  for (i in seq_along(bin_counts$y)) {
    cat("Bin:", i, ", Count:", bin_counts$count[i], "\n")
  }

  # Annotate which FOV the data belongs to
  myplot$fov = FOV
  # make myplot df unique for each fov ##########
  # Creating a new object with a dynamic name
  dfTitle <- "UMAP1_midline_data"
  dfName <- paste(dfTitle, FOV, sep = "_")
  assign(dfName, myplot)

  print("If spatial midline looks off, please adjust 'degs' parameter!")
  message('Please save the dataframe with FOV identifier, such as:  ',
          dfName)

  #change default assay back to whatever it was before running function
  SeuratObject::DefaultAssay(object) <- current_default_assay

  # Return the resulting df (for the FOV specified)
  return(get(dfName))
}

#' Plots UMAP_1 embeddings values for cells and their distance away from the computed spatial midline (supports multiple FOVs)
#' @author Margarita Kapustina
#'
#' @description Plot UMAP_1 embedding values for cells and their distance away from Spatial Midline across multiple FOVs.
#'
#' @param UMAP1MidlineData List of dataframes generated with getUMAP1_MidlineData() across multiple FOVs
#' @param binNumber Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth
#' @param binwidth Width of bins for gene expression averaging (lines). Suggested binwidth = 40 for 40micron bins
#' @param save_plot Option to save plot as .eps in working directory (TRUE, FALSE)
#'
#' @return Pooled dataframe containing Cell IDs, coordinates (X,Y), UMAP_1 embedding values,
#'    and computed distance away from spatial midline of each cell, across all FOVs in data
#'
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import tibble
#'
#' @examples
#' pooled_UMAP1_vsMidline = plotUMAP1_vsMidline(list(UMAP1_midline_data_fov, UMAP1_midline_data_X1fov))
#' @export

plotUMAP1_vsMidline <- function(UMAP1MidlineData,
                                binNumber = 7, binwidth =40,
                                save_plot = FALSE) {
  # Use do.call to rbind the objects in the list
  pooled_UMAP_df <- do.call(rbind, UMAP1MidlineData)

  # ggplot2 Visualization
  #show UMAP_1 embeddings as a line graph, as a distance away from midline
  p1 = ggplot2::ggplot(pooled_UMAP_df,ggplot2::aes(x = dLine)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count / 100)), bins = binNumber, boundary= 0, colour = 'green')+
    ggplot2::coord_cartesian(xlim = c(-500, 500), ylim =c(-6,8))+
    ggplot2::xlab('Distance away from Spatial Midline')+
    ggplot2::ylab('Cells per bin/100')+
    ggplot2::stat_summary_bin(ggplot2::aes(x = dLine, y= UMAP_1),fun = mean, binwidth =binwidth, geom = 'line', colour = 'blue') +
    #stat_bin(bins = 7, boundary = 0,geom = 'text', aes(label = ..count../100), vjust = 0, size = 3, color = 'black')+
    ggplot2::ggtitle('UMAP_1 embeddings across Spatial Midline')
  print(p1)

  # Calculate the number of cells in each bin
  bin_counts <- ggplot2::ggplot_build(p1)$data[[1]]
  cat("Number of cells in each bin (note: plotted are counts/100):
Beginning on left-most bin...\n")
  for (i in seq_along(bin_counts$y)) {
    cat("Bin:", i, ", Count:", bin_counts$count[i], "\n")
  }

  # saving
  if (save_plot) {
    print('Saving plot...')
    dev.copy2eps(file = 'UMAP1_embeddings_vs_midline.eps')
    print("Plot saved to working directory.")}

  return(pooled_UMAP_df)
}
