#' Generate Gene Expression values for cells and their distance away from Spatial Midline
#' @author Margarita Kapustina
#'
#' @description Generates gene expression data away from the computed spatial midline. Supports one FOV at a time.
#'
#' @param object a Xenium object
#' @param genes Character vector of gene names to fetch expression data for
#' @param FOV FOV to extract coordinates from
#' @param degs User-defined angle, defined in degrees, that intersects centrepoint in midline computation
#'
#' @return Dataframe containing Cell IDs, coordinates (X,Y), gene expression counts for specififed gene(s),
#'    and computed distance away from spatial midline of each cell, in specified FOV.
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import tibble
#'
#' @examples
#' geneExpression_X1fov <- getExpressionvsMidline(AD, FOV = "X1fov", genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'), degs = 45)
#' geneExpression_fov <- getExpressionvsMidline(AD, FOV = "fov", genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'), degs = 50)
#' @export

getExpressionvsMidline <- function(object, FOV = 'fov', genes, degs) {

  #get current default assay to then set back later
  current_default_assay = SeuratObject::DefaultAssay(object)

  # Set default FOV
  SeuratObject::DefaultFOV(object) <- FOV

  # Get tissue coordinates for the set default FOV only
  df <- Seurat::GetTissueCoordinates(object)

  # Fetch gene expression data
  exp <- Seurat::FetchData(object, genes, layer = 'counts')
  exp = tibble::rownames_to_column(exp, var = "cell")

  # Get expression values for only FOV of interest
  cells_to_get <- df[["cell"]]
  exp_subset <- exp[exp$cell %in% cells_to_get, ]

  # Make an expression and coordinate table
  expressionTable_fov <- cbind(exp_subset, df)
  expressionTable_fov <- expressionTable_fov[, !duplicated(colnames(expressionTable_fov))]

  # Change column names to capitals (works for all number of genes)
  names(expressionTable_fov)[names(expressionTable_fov) == 'x'] <- 'X'
  names(expressionTable_fov)[names(expressionTable_fov) == 'y'] <- 'Y'

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

  # Compute center of coordinates used for analysis (i.e., per FOV)
  centre = expressionTable_fov %>% .getCentre(doMean = FALSE)
  # Run distance operation, make a new data frame "myplot"
  rads <- degs * pi / 180
  m <- tan(rads)
  dLine <- (-m * centre$dX + centre$dY) / sqrt(1 + m^2)
  myplot <- cbind(centre, dLine)
  #plot out UMAP_1 embedding with nucleus midline
  gg <- ggplot2::ggplot(myplot,aes(x=dX,y=dY))
  gg <- gg + ggplot2::geom_point()
  gg <- gg + ggplot2::geom_abline(slope=m)
  gg <- gg + ggplot2::coord_fixed() +ggtitle('Midline of Cells')
  print(gg)

  # Annotate which FOV the data belongs to
  myplot$fov = FOV
  # make myplot df unique for each fov ##########
  # Creating a new object with a dynamic name
  dfTitle <- "geneExpression"
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


#' Plots gene expression data away from the computed spatial midline (supports multiple FOVs)
#' @author Margarita Kapustina
#'
#' @description Plot gene expression data for cells and their distance away from Spatial Midline across multiple FOVs.
#'
#' @param geneExpressionData List of dataframes generated with getExpressionvsMidline() across multiple FOVs
#' @param genes Character vector of gene names to their plot expression data
#' @param binNumber Number of bins to use for histogram (bars). Suggested binNumber = total # cells in pooled data/binwidth
#' @param binwidth Width of bins for gene expression averaging (lines) Suggested binwidth = 40 for 40micron bins
#' @param save_plot Option to save plot as .eps in working directory (TRUE, FALSE)
#'
#' @return Pooled dataframe containing Cell IDs, coordinates (X,Y), gene expression counts for specififed gene(s),
#'    and computed distance away from spatial midline of each cell, across all FOVs in data
#'
#' @import ggplot2
#' @import Seurat
#' @import tibble
#'
#' @examples
#' pooled_GeneExpression = plotGeneExpressionVsMidline(geneExpressionData = list(geneExpression_fov, geneExpression_fX1ov), genes = c("Epha1", "Cntnap4", 'Pcp4', 'Slc17a7'))
#' @export


plotGeneExpressionVsMidline <- function(geneExpressionData, genes,
                                        binNumber = 7, binwidth = 40,
                                        save_plot = FALSE) {
  # Use do.call to rbind the objects in the list
  pooled_df <- do.call(rbind, geneExpressionData)
  #create plot
  p1 <- ggplot2::ggplot(pooled_df, aes(x = dLine)) +
    geom_histogram(aes(y = stat(count / 100)), bins = binNumber, boundary = 0, colour = 'blue') +
    coord_cartesian(ylim = c(0, 5), xlim = c(-500, 500)) +
    xlab('Distance away from spatial midline') +
    ylab(paste('Bars: Number of cells per binned distance/100 \n Lines: Gene expression averaged per bin for genes:\n', paste(genes, collapse = ", "))) +
    ggtitle('Average gene expression values vs Spatial midline')

  # Define a color palette (adjust as needed)
  color_palette <- c("red", "green", "pink", "black", "orange")
  # Create a mapping from gene names to colors
  gene_color_mapping <- setNames(color_palette[1:length(genes)], genes)
  cat("Colour Scheme (Gene : colour):\n")
  # Add lines for each gene in genes with consistent colors
  for (gene in genes) {
    gene_color <- gene_color_mapping[[gene]]
    p1 <- p1 + ggplot2::stat_summary_bin(aes(x = dLine, y = .data[[gene]]), fun = mean, binwidth = binwidth, geom = 'line', colour = gene_color)
    # Print message indicating which gene and color are being used
    cat(gene, ":", gene_color, "\n")}
  print("Generating plot...")
  print(p1)

  # Calculate the number of cells in each bin
  bin_counts <- ggplot_build(p1)$data[[1]]
  cat("Number of cells in each bin (note: plotted are counts/100):
Beginning on left-most bin...\n")
  for (i in seq_along(bin_counts$y)) {
    cat("Bin:", i, ", Count:", bin_counts$count[i], "\n")
  }
  # saving
  if (save_plot) {
    print('Saving plot...')
    dev.copy2eps(file = 'BoxPlots_SCTcounts.eps')
    print("Plot saved to working directory.")}

  return(pooled_df)
}
