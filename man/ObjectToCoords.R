#' Get cell coordinates from a Xenium object
#' @author Margarita Kapustina
#'
#' @description Create a dataframe with cell ID coordinates and cluster IDs from a Xenium Seurat object.
#'
#' @param object Xenium object
#' @param thisFOV Name of your FOV to extract coordinates from
#' @param threshold_y Boundary that you define to compute cells that are above and below this boundary
#' @param angle_adjust Option to adjust the angle of your cell coordinates (TRUE, FALSE)
#' @param thetda_deg Specify the degrees you wish to rotate your slice by
#' @param flip_x_coordinates Option to flip your slice in the horizontal plane (TRUE, FALSE)
#'
#' @return A dataframe with cell ID, coordinates, cluster ID, and a plot with UMAP_1 embedding values and the boundary that you define
#'
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import cowplot
#'
#' @examples
#' distanceData_X2fov = object_FOV_to_coordinates(object = a, thisFOV = 'X2fov', threshold_y = 2694, angle_adjust = FALSE, theta_deg = 0, flip_x_coordinates = FALSE)
#' @export
object_FOV_to_coordinates <- function(object = object, thisFOV = 'X9fov', threshold_y = 2000, angle_adjust = FALSE, theta_deg = 0, flip_x_coordinates = FALSE) {
  #set default fov
  ####################################
  DefaultFOV(object) = thisFOV
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
  # Get tissue coordinates (for slice Xfov3)
  df <- Seurat::GetTissueCoordinates(object)
  # Merge tissue coordinates (FOV of interest only) with the dataframe containing 1D UMAP embeddings
  cells_to_get <- df[["cell"]]
  cell.embeddings.df_subset <- cell.embeddings.df[cell.embeddings.df$Cell_ID %in% cells_to_get, ]
  names(cell.embeddings.df_subset)[1] <- 'cell'
  coordsEmbeddingsdf_oneFOV <- merge(df, cell.embeddings.df_subset, by = 'cell')
  # Angle of rotation in degrees if angle_adjust is TRUE
  if (angle_adjust) {
    theta_rad <- theta_deg * (pi / 180) #deg to radians
    newX <- coordsEmbeddingsdf_oneFOV$x * cos(theta_rad) + coordsEmbeddingsdf_oneFOV$y * sin(theta_rad)
    newY <- -coordsEmbeddingsdf_oneFOV$x * sin(theta_rad) + coordsEmbeddingsdf_oneFOV$y * cos(theta_rad)
    coordsEmbeddingsdf_oneFOV$x <- newX
    coordsEmbeddingsdf_oneFOV$y = newY
  }

  if (flip_x_coordinates) {
    # Calculate the maximum x-coordinate value
    max_x <- max(coordsEmbeddingsdf_oneFOV$x)
    coordsEmbeddingsdf_oneFOV$x <- max_x - coordsEmbeddingsdf_oneFOV$x # Flip the x-coordinates
  }
  ##############################
  # Define the threshold x-value for each indv. fov ##############change for each FOV!
  #add what the threshold was
  coordsEmbeddingsdf_oneFOV$threshold <- threshold_y
  # Create a new column based on the condition
  coordsEmbeddingsdf_oneFOV$above_or_below <- ifelse(coordsEmbeddingsdf_oneFOV$y > threshold_y, "Above", "Below")
  #add what the fov was
  coordsEmbeddingsdf_oneFOV$fov = thisFOV

  message('Your threshold is...')
  print(threshold_y)

  message('Plotting your final coordinates, UMAP_1 embeddings and threshold...')
  #plot out UMAP_1 embedding with nucleus midline####################
  gg <- ggplot(coordsEmbeddingsdf_oneFOV,aes(x=x,y=y))
  gg <- gg + geom_point(aes(colour= UMAP_1), size = 1.5)
  gg <- gg + geom_hline(yintercept = threshold_y, linetype = "dashed", color = "red")
  gg <- gg + coord_fixed() +ggtitle(thisFOV,threshold_y)+scale_color_viridis_c()
  print(gg) ##############print  for each FOV!
  # Construct the file name dynamically based on the value of this
  message('Plot will be saved as:')
  file_name <- paste(thisFOV, ".eps", sep = "")
  print(file_name)
  ggsave(file_name, gg)
  # Extract cluster IDs from Seurat object ###########
  cluster_ids <- as.data.frame(Idents(object))
  colnames(cluster_ids) <- c("cluster_ids")  # Renaming column for better readability
  cluster_ids$cell <- rownames(cluster_ids)# Adding cell ID column based on row names
  rownames(cluster_ids) = NULL
  #View(cluster_ids)
  # Merge cluster IDs into the existing dataframe based on cell IDs############
  merged_data <- merge(coordsEmbeddingsdf_oneFOV, cluster_ids, by = "cell", all.x = TRUE)
  #name dataframe with unique identifier, and eename the object dynamically
  uniqueFOV_name <- paste("distanceData_", thisFOV, sep = "")
  assign(uniqueFOV_name, merged_data)
  ############
  message('Name this file as:')
  print(uniqueFOV_name)
  return(merged_data)
  ################################# END ANALYSIS HERE
}
