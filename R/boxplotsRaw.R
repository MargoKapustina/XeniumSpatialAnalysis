#' Xenium Gene Box Plots - Raw Counts
#' @author Margarita Kapustina
#'
#' @description Create a boxplot denoting expression of gene(s) of interest across all clusters. Xenium assay counts (raw counts) are used.
#'
#' @param object Xenium object
#' @param genes Character vector of gene names to plot
#' @param ncol Specify number of columns in plot
#' @param save_plot Option to save plot as .eps in working directory (TRUE, FALSE)
#'
#' @return A grid of boxplots showing raw gene expression across all clusters
#'
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import cowplot
#'
#' @examples
#' #create vector ie my_genes= c('C1ql2','Slc17a7', 'Gng13')
#' BoxPlots = XeniumBoxPlotRaw(object = atn, genes = my_genes)
#' @export

XeniumBoxPlotRaw = function(object, genes, ncol = 4, save_plot = FALSE) {

  #get current default assay to then set back later
  current_default_assay = SeuratObject::DefaultAssay(object)
  #change default assay to Xenium, otherwise default is SCT
  SeuratObject::DefaultAssay(object) <- "Xenium"

  suppressWarnings({ #supresses following warning:
    # Computation failed in `stat_ydensity()`
    #Caused by error in `density.default()`:
    # ! 'bw' is not positive.

    # Function to generate a customized violin plot
    Plot = function(theGenes) {
      Seurat::VlnPlot(object, features = theGenes, pt.size = 0, layer = "counts", adjust = 0) +
        ggplot2::theme(legend.position = "nolegend") +
        ggplot2::geom_boxplot(outlier.size = 0, outlier.stroke = 0, lwd = 0.3) +
        ggplot2::theme(axis.line = element_line(linewidth = 0.5), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6),
              axis.title = element_blank(), title = element_text(size = 6),
              plot.background = element_blank())
    }
    print('Generating boxplots...')

    # Initialize an empty list to store individual plots
    plots <- list()

    # Loop through the gene vector and call the Plot function for each gene
    for (i in seq_along(genes)) {
      plotName <- paste0("p", i)
      assign(plotName, Plot(theGenes = genes[i]))

      # Add the plot to the list
      plots[[i]] <- get(plotName)
    }

    # Create a grid of plots using plot_grid
    grid <- cowplot::plot_grid(plotlist = plots, ncol = ncol)

    # Print the combined plot
    print(grid)

    # saving
    if (save_plot) {
      print('Saving plot...')
      dev.copy2eps(file = 'BoxPlots_rawCounts.eps')
      print("Plot saved to working directory.")}

    #change default assay back to whatever it was before running function
    SeuratObject::DefaultAssay(object) <- current_default_assay

    # Return the combined plot
    return(grid)

    # Code that may produce warnings
  })
}
