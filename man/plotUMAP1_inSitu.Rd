% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotUMAP1_inSItu.R
\name{plotUMAP1_inSitu}
\alias{plotUMAP1_inSitu}
\title{Xenium 1-D UMAP Embeddings Plotted in Situ}
\usage{
plotUMAP1_inSitu(
  object,
  FOV = "X1fov",
  EmbeddingsPlotTitle = "UMAP_1 embedding for all\\n cells in subset",
  HistogramPlotTitle = "Distribution of UMAP_1 embedding\\n values all cells in subset",
  inSituPlotTitle = "UMAP_1 embedding in situ\\n for representative FOV",
  save_plot = FALSE
)
}
\arguments{
\item{object}{a Xenium object}

\item{FOV}{FOV to plot UMAP_1 embedding values in}

\item{EmbeddingsPlotTitle}{title for UMAP_1 Embeddings plot}

\item{HistogramPlotTitle}{title for UMAP_1 Embeddings Histogram plot}

\item{inSituPlotTitle}{title for UMAP_1 Embeddings in situ plot}

\item{save_plot}{Option to save plot as .eps in working directory (TRUE, FALSE)}
}
\value{
Dataframe containing Cell IDs, coordinates (X,Y) and computed UMAP_1 embeddings values for cells in specified FOV.
#Plots 1D UMAP, corresponding histogram of UMAP_1 embeddings values, and UMAP_1 embeddings values plotted _in situ_.
}
\description{
Plots: 1-dimensional UMAP, the corresponding histogram of UMAP_1 embedding values, and
   the UMAP_1 embedding values in situ within specified FOV.
}
\examples{
AD_df = plotUMAP1_inSitu(object = AD, FOV = 'fov')
}
\author{
Margarita Kapustina
}
