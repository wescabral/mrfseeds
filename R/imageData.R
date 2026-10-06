#' Image Data Class
#'
#' An R6 class that handles image data with support for both discrete and continuous values.
#' Provides methods for visualization and matrix manipulation.
#'
#' @details
#' The `imageData` class wraps a numeric matrix and provides utilities for:
#' - Automatic detection of discrete vs continuous data
#' - Visualization via ggplot2
#' - Management of image dimensions and value ranges
#'
#' @examples
#' \dontrun{
#' # Create an imageData object from a discrete matrix
#' z <- matrix(sample(0:2, 100, replace = TRUE), nrow = 10, ncol = 10)
#' img <- imageData$new(z)
#'
#' # Plot the image
#' img$plot()
#'
#' # Create a continuous image
#' y <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' img_cont <- imageData$new(y, continuous = TRUE)
#' img_cont$plot()
#'
#' # Update matrix - noisy example
#' new_z <- img$matrix + rnorm(prod(dim(img$matrix)), sd = 0.1)
#' img$setMatrix(new_z)
#' img$plot()
#' }
#'
#' @export
imageData <- R6::R6Class("imageData",
  public = list(

    #' @description
    #' Initialize an imageData object
    #'
    #' @param dataMatrix A numeric matrix representing image data
    #' @param continuous Logical. If `NULL`, automatically detects if data is continuous
    #'   or discrete. If `TRUE`, treats all values as continuous. If `FALSE`, treats
    #'   as discrete/categorical.
    #'
    #' @return An imageData object
    initialize = function(dataMatrix, continuous = NULL) {
      if(!is.matrix(dataMatrix)){
        stop("dataMatrix must be a matrix.")
      }
      private$.dataMatrix <- dataMatrix
      private$.dim <- dim(dataMatrix)

      if(is.null(continuous)){
        if(is.numeric(dataMatrix)){
          vals <- as.vector(dataMatrix)
          epsilon <- .Machine$double.eps^0.5
          if(all(abs(vals - round(vals)) < epsilon, na.rm = TRUE)){
            private$.continuous <- FALSE
          } else {
            private$.continuous <- TRUE
          }
        } else {
          private$.continuous <- FALSE
        }
      } else {
        private$.continuous <- isTRUE(continuous)
      }

      if(private$.continuous){
        private$.valueRange <- range(as.vector(dataMatrix), na.rm = TRUE)
        private$.valueSet <- NULL
        private$.ncolors <- NA
      } else {
        private$.valueSet <- unique(as.vector(dataMatrix))
        private$.ncolors <- length(private$.valueSet)
        private$.valueRange <- range(as.vector(dataMatrix), na.rm = TRUE)
      }
    },

    #' @description
    #' Update the image matrix with new data
    #'
    #' @param newdataMatrix A numeric matrix or vector. If a vector, it will be
    #'   reshaped to match the dimensions of the current matrix.
    #' @param continuous Logical. If `NULL`, automatically detects if new data is
    #'   continuous or discrete. Otherwise, overrides automatic detection.
    #'
    #' @return Invisibly returns `self` for method chaining
    setMatrix = function(newdataMatrix, continuous = NULL) {
      if(is.vector(newdataMatrix) && length(newdataMatrix) == prod(private$.dim)){
        newdataMatrix <- matrix(newdataMatrix, nrow = private$.dim[1], ncol = private$.dim[2])
      }
      if(!is.matrix(newdataMatrix)){
        stop("dataMatrix must be a matrix or a vector with compatible length as the prior matrix.")
      }
      private$.dataMatrix <- newdataMatrix
      private$.dim <- dim(newdataMatrix)

      if(is.null(continuous)){
        if(is.numeric(private$.dataMatrix)){
          vals <- as.vector(private$.dataMatrix)
          epsilon <- .Machine$double.eps^0.5
          if(all(abs(vals - round(vals)) < epsilon, na.rm = TRUE)){
            private$.continuous <- FALSE
          } else {
            private$.continuous <- TRUE
          }
        } else {
          private$.continuous <- FALSE
        }
      } else {
        private$.continuous <- isTRUE(continuous)
      }

      if(private$.continuous){
        private$.valueRange <- range(as.vector(private$.dataMatrix), na.rm = TRUE)
        private$.valueSet <- NULL
        private$.ncolors <- NA
      } else {
        private$.valueSet <- unique(as.vector(private$.dataMatrix))
        private$.ncolors <- length(private$.valueSet)
        private$.valueRange <- range(as.vector(private$.dataMatrix), na.rm = TRUE)
      }
      invisible(self)
    },

    #' @description
    #' Generate a ggplot2 visualization of the image
    #'
    #' @details
    #' For discrete data, uses the "Set3" color palette with categorical fill.
    #' For continuous data, uses a gradient from white to steelblue.
    #'
    #' @return A ggplot2 object
    plot = function(){
      df <- data.frame(x = as.vector(row(private$.dataMatrix)),
                       y = as.vector(col(private$.dataMatrix)),
                       value = as.vector(private$.dataMatrix))

      base_plot <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y))

      if(private$.continuous){
        base_plot <- base_plot +
          ggplot2::geom_tile(ggplot2::aes(fill = value)) +
          ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
          ggplot2::labs(fill = NULL)
      } else {
        base_plot <- base_plot +
          ggplot2::geom_tile(ggplot2::aes(fill = as.factor(value))) +
          ggplot2::scale_fill_brewer(palette = "Set3") +
          ggplot2::labs(fill = NULL)
      }

      base_plot +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 5)))
    }
  ),

  active = list(

    #' @description
    #' Get the underlying matrix
    #'
    #' @return The numeric matrix
    matrix = function() {
      private$.dataMatrix
    },

    #' @description
    #' Get the dimensions of the image
    #'
    #' @return A numeric vector of length 2: c(nrow, ncol)
    dim = function() {
      private$.dim
    },

    #' @description
    #' Get the set of unique values (for discrete data)
    #'
    #' @return A vector of unique values, or `NULL` for continuous data
    valueSet = function() {
      private$.valueSet
    },

    #' @description
    #' Get the number of colors/classes (for discrete data)
    #'
    #' @return An integer representing the number of unique values, or `NA` for continuous data
    ncolors = function() {
      private$.ncolors
    },

    #' @description
    #' Check if the data is continuous
    #'
    #' @return Logical
    continuous = function() {
      private$.continuous
    },

    #' @description
    #' Get the range of values in the image
    #'
    #' @return A numeric vector of length 2: c(min, max)
    valueRange = function() {
      private$.valueRange
    }
  ),

  private = list(
    .dataMatrix = NULL,
    .dim = NULL,
    .valueSet = NULL,
    .ncolors = NULL,
    .continuous = FALSE,
    .valueRange = NULL
  )
)
