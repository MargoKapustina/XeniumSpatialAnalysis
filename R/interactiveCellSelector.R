#' Interactive cell selector: lasso tool
#'
#' Launches a Shiny app with a spatial cell coordinate plot for interactively selecting cells within a region of interest in a given FOV.
#' Users can drag-to-select cells and save the selected cell IDs and corresponding filtered coordinates to the global environment.
#'
#' @param fov_label A character string specifying the field of view (FOV) label.
#' This is used to name the output variables saved to the global environment. 
#'
#' @param cell_coords_df A dataframe containing cell coordinate data, generated using \code{object_FOV_to_coordinates()}. 
#' The dataframe must include the following columns:
#' \itemize{
#'   \item \code{x} : x-coordinates of cells
#'   \item \code{y} : y-coordinates of cells
#'   \item \code{cluster_ids} : cluster identity for coloring
#'   \item \code{cell} : unique cell IDs
#' }
#' 
#' @param plot_height Character string specifying plot height. 
#' This controls the vertical display size of the interactive FOV visualization. 
#' Must include CSS units (e.g. "750px"). The default is 750px.
#' 
#' @param pt.size Numeric value controlling point size in the scatter plot (default 0.8)
#'
#' @details
#' The Shiny app displays an interactive spatial plot of cells colored by cluster identity.
#' Drag to select cells of interest. Upon selection, click the "Save Selected Cell IDs"
#' button, which will save two objects in the global environment:
#' \itemize{
#'   \item \code{cellid_<your_fov_label>} : a vector of selected cell IDs (e.g., distanceData_X1fov). Use this to subset cells via: \code{subsetObj <- subset(myObj, cells = cellid_<your_fov_label>)} 
#'   \item \code{<your_fov_label>} : a filtered dataframe containing only the selected cells. This object will be saved using your specified \code{fov_label} as its name.
#'   (e.g., x1fov). 
#'}
#' Once you have clicked the "Save Selected Cell IDs" button, close the window. If no cells are selected, a notification will be displayed. 
#' The saved cell ID vectors can be combined across FOVs to subset an entire region of interest using: 
#' \code{combinedCellids <- c(distanceData_X1fov, distanceData_X2fov)
#' SubsetObj <- subset(object, cells = combinedCellids)}
#' 
#' @return Launches an interactive Shiny application. No value is returned.
#'
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#' @import shiny
#' @importFrom plotly ggplotly layout event_data
#'
#' @examples
#' interactiveCellSelector("x1fov", distanceData_X1fov)
#' @export

interactiveCellSelector <- function(
    fov_label,
    cell_coords_df,
    plot_height = "750px",
    pt.size = 0.8
) {
  
  # update to df 
  df = cell_coords_df
  # check data quality
  required_cols <- c("x", "y", "cluster_ids", "cell")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ".\nYou can generate df with cell coordinates using object_FOV_to_coordinates()"
    )
  }
  
  # strict px validation
  if (!grepl("^[0-9]+(px)$", plot_height)) {
    stop("plot_height must be a character string with px units. Example: '750px'.")
  }
  
  # plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = x,
      y = y,
      color = as.factor(cluster_ids),
      key = cell
    )
  ) +
    ggplot2::geom_point(size = pt.size) +
    ggplot2::coord_fixed() +
    ggplot2::labs(color = "Cluster")
  
  pp <- plotly::ggplotly(p) |>
    plotly::layout(dragmode = "lasso")
  
  # UI
  ui <- shiny::fluidPage(
    
    shiny::tags$head(
      shiny::tags$style(
        shiny::HTML("
          .app-title {
            font-size: 18px;
            font-weight: 600;
            margin-bottom: 10px;
          }
        ")
      )
    ),
    
    shiny::div(
      class = "app-title",
      "Interactive Cell Selector"
    ),
    
    plotly::plotlyOutput("plt", height = plot_height, width = "100%"),
    shiny::verbatimTextOutput("sel"),
    
    shiny::actionButton(
      "save_btn",
      "Save Selected Cell IDs",
      class = "btn btn-primary"
    )
  )
  
  # server
  server <- function(input, output, session) {
    
    output$plt <- plotly::renderPlotly(pp)
    
    selected_ids <- shiny::reactive({
      s <- plotly::event_data("plotly_selected")
      if (is.null(s)) return(NULL)
      s$key
    })
    
    output$sel <- shiny::renderPrint({
      ids <- selected_ids()
      
      if (is.null(ids)) {
        cat("No cells selected yet")
      } else {
        cat("Preview of cell IDs:\n")
        utils::tail(ids)
      }
    })
    
    shiny::observeEvent(input$save_btn, {
      
      ids <- selected_ids()
      
      if (is.null(ids)) {
        
        shiny::showNotification(
          "No selection!",
          type = "error"
        )
        
      } else {
        
        id_var <- paste0("cellid_", fov_label)
        df_var <- fov_label
        
        # save IDs
        base::assign(id_var, ids, envir = .GlobalEnv)
        
        # filter df
        filtered_df <- df[df$cell %in% ids, , drop = FALSE]
        base::assign(df_var, filtered_df, envir = .GlobalEnv)
        
        message(
          "Saved: ", id_var, " and ", df_var,
          " | ", length(ids), " cells selected"
        )
        
        shiny::showNotification(
          shiny::HTML(
            paste0(
              "Saved: ", id_var, " and ", df_var, "<br>",
              "[Selection contains: ", length(ids), " cells]"
            )
          ),
          type = "message"
        )
      }
    })
  }
  
  shiny::shinyApp(ui, server)
}