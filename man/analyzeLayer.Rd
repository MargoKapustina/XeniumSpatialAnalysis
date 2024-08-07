% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyzeLayer.R
\name{analyzeLayer}
\alias{analyzeLayer}
\title{Calculate the normalized distance along a cortical layer and away from a boundary for cell coordinates}
\usage{
analyzeLayer(thisFOV)
}
\arguments{
\item{thisFOV}{Name of your FOV to define your boundary within}
}
\value{
A dataframe with cell ID, coordinates, cluster ID, normalized distance along a cortical layer, and distance away from your defined boundary.
}
\description{
Compute distance along a cortical layer, and the distance away from a defined boundary, using distance data computed with object_FOV_to_coordinates function. Note: please save your data computed object_FOV_to_coordinates named as 'distance_data_all'. If you have multiple FOVs run through the object_FOV_to_coordinates function, merge them and rename the object. If you have only a single FOV, simply rename the object.
}
\examples{
#before running the function, please save merged FOV-specific data computed with the object_FOV_to_coordinates function, as distance_data_all
#example: distance_data_all = rbind(distanceData_X2fov, distanceData_X3fov)
#Please noe that distance_data_all needs to be within your environment to run analyzeLayer
cellDataX3fov <- analyzeLayer("X3fov")
}
\author{
Margarita Kapustina
}
