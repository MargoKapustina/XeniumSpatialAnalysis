#' Calculate the normalized distance along a cortical layer and away from a boundary for cell coordinates
#' @author Margarita Kapustina
#'
#' @description Compute distance along a cortical layer, and the distance away from a defined boundary, using distance data computed with object_FOV_to_coordinates function. Note: please save your data computed object_FOV_to_coordinates named as 'distance_data_all'. If you have multiple FOVs run through the object_FOV_to_coordinates function, merge them and rename the object. If you have only a single FOV, simply rename the object.
#'
#' @param thisFOV Name of your FOV to define your boundary within
#'
#' @return A dataframe with cell ID, coordinates, cluster ID, normalized distance along a cortical layer, and distance away from your defined boundary.
#'
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import cowplot
#'
#' @examples
#' #before running the function, please save merged FOV-specific data computed with the object_FOV_to_coordinates function, as distance_data_all
#' #example: distance_data_all = rbind(distanceData_X2fov, distanceData_X3fov)
#' #Please noe that distance_data_all needs to be within your environment to run analyzeLayer
#' cellDataX3fov <- analyzeLayer("X3fov")
#' @export

#run on distance data all (get from objectToCoords function)
analyzeLayer<- function(thisFOV) {
  #select FOV to get subset of data
  cell_data <- distance_data_all[distance_data_all$fov == thisFOV, ]

  # Plot out the cells in space
  # ggplot(cell_data, aes(x = x, y = y)) +
  #  geom_point()+
  #  coord_cartesian(ylim = c(0, 5000), xlim = c(0,5000))+
  #  coord_fixed()  # Set y-axis limits

  #doesnt mattter if aspect ratio is off, continue as normal the plot will re-adjust later when using ggplot2
  par(pty="s")
  plot(cell_data$x, cell_data$y, type = "p", xlab = "x", ylab = "y")

  message('Please click around 25-30 points to set major landmarks in plot screen...')
  message('When done hit ESC on keyboard.')

  # Prompt user to click points on the plot
  points <- locator() #hit ESC when done on keyboard

  # Extract x and y coordinates of clicked points
  x_coords <- points$x
  y_coords <- points$y
  points = as.data.frame(points)

  # ggplot(points, aes(x = x, y = y)) +
  #  geom_point(color = "blue")+
  # coord_cartesian(ylim = c(0, 5000), xlim = c(0,5000)) +
  #geom_point(data = cell_data,  aes(x = x, y = y), color = "red") + coord_fixed()

  # Plot the points and create the path
  #  ggplot(points, aes(x = x, y = y)) +
  #   geom_path(color = "blue") +
  #  geom_point(color = "red")


  ######### analysis v1
  # Function to calculate the spline between neighboring points
  calculate_spline <- function(x1, y1, x2, y2, n_points = 100) {
    # Generate x values for the spline
    x_values <- seq(x1, x2, length.out = n_points)

    # Interpolate y values using a spline
    y_values <- spline(x = c(x1, x2), y = c(y1, y2), xout = x_values)$y

    return(data.frame(x = x_values, y = y_values))
  }

  # Initialize a list to store spline coordinates
  spline_segments <- list()

  # Iterate through points and calculate spline segments
  for (i in 1:(nrow(points) - 1)) {
    # Calculate spline segment between neighboring points
    segment <- calculate_spline(points$x[i], points$y[i], points$x[i + 1], points$y[i + 1])
    spline_segments[[i]] <- segment
  }

  # Combine spline segments into a single DataFrame
  spline_data <- do.call(rbind, spline_segments)

  # Plot the spline data
  # plot(spline_data$x, spline_data$y, type = "p",
  #     col = "blue", xlab = "X", ylab = "Y", main = "Spline Through Nearest Points")

  project_point_onto_segment = function(px, py, x1, y1, x2, y2) {
    dx <- x2 - x1
    dy <- y2 - y1
    t <- ((px - x1) * dx + (py - y1) * dy) / (dx^2 + dy^2)
    t <- pmax(0, pmin(1, t))
    projected_x <- x1 + t * dx
    projected_y <- y1 + t * dy
    return(c(projected_x, projected_y))
  }
  # Initialize a vector to store the distances
  distances <- numeric(nrow(cell_data))

  # Plot the cell_data points, spline_data, and distances
  par(pty="s")
  plot(cell_data$x, cell_data$y, type = "p", col = "red", pch = 16, xlab = "X", ylab = "Y", main = "Cell coordinates, estimated L6b boundary and Distance to boundary")
  points(spline_data$x, spline_data$y, col = "blue")
  for (i in 1:nrow(cell_data)) {
    nearest_index <- which.min((spline_data$x - cell_data$x[i])^2 + (spline_data$y - cell_data$y[i])^2)
    if (nearest_index == 1) {
      projected <- project_point_onto_segment(cell_data$x[i], cell_data$y[i],
                                              spline_data$x[1], spline_data$y[1],
                                              spline_data$x[2], spline_data$y[2])
    } else if (nearest_index == nrow(spline_data)) {
      projected <- project_point_onto_segment(cell_data$x[i], cell_data$y[i],
                                              spline_data$x[nrow(spline_data) - 1], spline_data$y[nrow(spline_data) - 1],
                                              spline_data$x[nrow(spline_data)], spline_data$y[nrow(spline_data)])
    } else {
      projected <- project_point_onto_segment(cell_data$x[i], cell_data$y[i],
                                              spline_data$x[nearest_index - 1], spline_data$y[nearest_index - 1],
                                              spline_data$x[nearest_index], spline_data$y[nearest_index])
    }
    segments(projected[1], projected[2], cell_data$x[i], cell_data$y[i], col = "green")  # Reversed the order of endpoints

    # Calculate the distance between the point and its projection onto the spline
    distances[i] <- sqrt((cell_data$x[i] - projected[1])^2 + (cell_data$y[i] - projected[2])^2)
  }

  # Add the distances as a new column in the cell_data dataframe
  cell_data$distance_to_spline <- distances

  # Add a legend
  legend("bottomright", legend = c("Cell Data", "Spline Data", "Distance to Spline"), col = c("red", "blue", "green"), pch = c(16, NA, NA), lty = c(NA, 1, 1))
  #plot: colour by distance away
  # Define the color palette using colorRampPalette
  my_palette <- colorRampPalette(c("magenta", "green"))(nrow(cell_data))

  # plot(cell_data$x, cell_data$y,
  #     col = my_palette[cell_data$distance_to_spline], pch = 16, xlab = "X", ylab = "Y",
  #    main = "Color-coded Cells by Distance to Spline")

  ##normalize the distance to spline stored in cell_data$distance_to_spline
  cell_data$normalized_distance_to_spline <- (cell_data$distance_to_spline - min(cell_data$distance_to_spline)) / (max(cell_data$distance_to_spline) - min(cell_data$distance_to_spline))

  ##ggplot version
  title <- paste0("Color-coded Cells by Normalized Distance to Spline: ", thisFOV)
  p2 =  ggplot(cell_data, aes(x = x, y = y)) +
    geom_point(aes(color = normalized_distance_to_spline), pch = 16) +
    scale_color_viridis_c() +
    labs(x = "X", y = "Y", title = title) +
    coord_fixed()
  # Define the file path and name for saving dynamically
  filename <- paste0("distancetoSpline_", thisFOV, ".eps")
  ggsave(filename = filename, plot = p2)
  ##or## as factor
  # plot(cell_data$x, cell_data$y,
  # col = my_palette[factor(cell_data$distance_to_spline)], pch = 16, xlab = "X", ylab = "Y",
  # main = "Color-coded Cells by Distance to Spline")

  ############### end here for now

  # Initialize vector to store distances
  cell_data$distance_along_spline <- numeric(nrow(cell_data))

  # Iterate over each point in cell_data
  for (i in 1:nrow(cell_data)) {
    # Calculate distances to all points in spline_data
    distances <- sqrt((spline_data$x - cell_data$x[i])^2 + (spline_data$y - cell_data$y[i])^2)

    # Find index of nearest point
    nearest_index <- which.min(distances)

    # Calculate distance along spline to nearest point
    if (nearest_index == 1) {
      # If nearest point is the first point in spline_data, set distance to 0
      cell_data$distance_along_spline[i] <- 0
    } else {
      # Calculate distance as sum of distances between consecutive points along the spline up to the nearest point
      cell_data$distance_along_spline[i] <- sum(sqrt(diff(spline_data$x[1:nearest_index])^2 + diff(spline_data$y[1:nearest_index])^2))
    }
  }

  ##normalize the distance to spline stored in cell_data$distance_along_spline
  cell_data$normalized_distance_along_spline <- (cell_data$distance_along_spline - min(cell_data$distance_along_spline)) / (max(cell_data$distance_along_spline) - min(cell_data$distance_along_spline))

  # Plot histogram of distances along the spline for all cells
  # Create the ggplot histogram
  title = paste0("Histogram of distances along the spline: ", thisFOV)
  p_hist <- ggplot(cell_data, aes(x = normalized_distance_along_spline)) +
    geom_histogram(color = "grey", fill = "violet", bins = 30) +
    labs(x = "Distance along spline", y = "Frequency", title = title)
  filename <- paste0("Histogram_distanceAlongSpline_", thisFOV, ".eps")
  ggsave(filename = filename, plot = p_hist)

  # Create the ggplot histogram
  title = paste0("Histogram of distances to the spline: ", thisFOV)
  p_hist2 <- ggplot(cell_data, aes(x = normalized_distance_to_spline)) +
    geom_histogram(color = "grey", fill = "gold", bins = 30) +
    labs(x = "Distance to spline", y = "Frequency", title = title)
  filename <- paste0("Histogram_distanceToSpline_", thisFOV, ".eps")
  ggsave(filename = filename, plot = p_hist2)
  ###############################

  # Color-code the cells by their distance along the spline
  my_palette <- colorRampPalette(c("red", "blue"))(nrow(cell_data))
  title <- paste0("Color-coded Cells by Normalized Distance Along Spline: ", thisFOV)
  p3 =  ggplot(cell_data, aes(x = x, y = y)) +
    geom_point(aes(color = normalized_distance_along_spline), pch = 16) +
    scale_color_gradientn(colors = my_palette) +
    labs(x = "X", y = "Y", title = title)+ coord_fixed()
  filename <- paste0("distanceAlongSpline_", thisFOV, ".eps")
  ggsave(filename = filename, plot = p3)

  #name dataframe with unique identifier, and name the object dynamically
  uniqueFOV_name <- paste("cellData", thisFOV, sep = "")
  assign(uniqueFOV_name, cell_data)
  ############
  message('Name this file as:')
  print(uniqueFOV_name)
  # Return the modified cell_data dataframe
  return(cell_data)
}
