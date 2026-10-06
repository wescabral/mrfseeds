#' Seed Model Class
#'
#' An R6 class that represents a Markov Random Field (MRF) model with centrality effects (seeds).
#' Generates synthetic images by sampling from the model using a Gibbs sampler.
#'
#' @details
#' The `seedModel` class implements an MRF model where:
#' - `alpha` is the uniformity parameter
#' - `beta` is the background class (color) penalty/favor
#' - `seeds` define one or more pixels with influence `delta` related to a certain `color`
#' - `ncolors` is the number of distinct color categories
#'
#' Images are generated using a Gibbs sampling algorithm that incorporates
#' the seeds as landmarks influencing the spatial pattern.
#'
#' @examples
#' \dontrun{
#' # Define seeds with positions and properties
#' real_seeds <- list(
#'   list(position = c(50, 50), color = 1, delta = 20),
#'   list(position = c(100, 100), color = 1, delta = 20)
#' )
#'
#' # Create model
#' mod <- seedModel$new(alpha = 0, beta = 1, seeds = real_seeds, ncolors = 2)
#'
#' # Generate an image
#' img <- mod$sampleImage(dim = c(150, 150), steps = 80)
#' img$plot()
#' }
#'
#' @export
seedModel <- R6::R6Class(
  "seedModel",
  public = list(

    #' @description
    #' Initialize a seed model
    #'
    #' @param alpha Numeric. The uniformity parameter. It controls the attraction between 
    #'   neighboring pixels. Higher parameters make it to smoother images.
    #' @param beta Numeric. The global background class/color strength. Higher values lead to
    #'   a image that favors the background color. 
    #' @param seeds A list of seed specifications. Each seed is a list with:
    #'   - `position`: numeric vector of length 2, c(x, y) coordinates
    #'   - `color`: integer, the preferred color/category at this seed
    #'   - `delta`: numeric, the influence radius of the seed (defines centrality effects). It has
    #'                to at least doubles `beta` for the influence radius visibility.
    #' @param ncolors Integer. The number of distinct color categories in the model.
    #'
    #' @return A seedModel object
    initialize = function(alpha, beta, seeds, ncolors){
      private$.alpha <- alpha
      private$.beta <- beta
      private$.seeds <- seeds
      private$.ncolors <- ncolors
    },

    #' @description
    #' Generate a synthetic image by sampling from the MRF model
    #'
    #' @details
    #' Uses a Gibbs sampler to generate images that respect the parameter given to the model.
    #'
    #' The number of sampling steps determines the mixing quality. More steps
    #' produce better samples but take longer.
    #'
    #' @param dim Numeric vector of length 2, c(nrow, ncol). Dimensions of the
    #'   generated image. Default: c(150, 150).
    #' @param steps Integer. Number of Gibbs sampling sweeps to perform.
    #'   Higher values produce better quality samples. Default: 60.
    #'
    #' @return An `imageData` object containing the sampled image
    sampleImage = function(dim = c(150, 150), steps = 60){
      z <- matrix(as.integer(sample(0L:(private$.ncolors - 1L), replace = TRUE, size = prod(dim))),
                  nrow = dim[1], ncol = dim[2])

      # Preparar matriz de seeds (n_seeds x 2)
      seedMatrix <- do.call(rbind, lapply(private$.seeds, function(s) s$position))

      # Preparar vetores de cores e deltas das seeds
      seedColors <- as.integer(sapply(private$.seeds, function(s) s$color))
      seedDeltas <- as.numeric(sapply(private$.seeds, function(s) s$delta))

      z <- gibbsSampler(z,
                        alpha = private$.alpha,
                        beta = private$.beta,
                        ncolors = private$.ncolors,
                        seeds = seedMatrix,
                        seedColors = seedColors,
                        seedDeltas = seedDeltas,
                        steps = steps)

      imageData$new(z)
    }

  ),

  active = list(

    #' @description
    #' Get the alpha parameter
    #'
    #' @return Numeric. The global uniformity parameter.
    alpha = function(){
      private$.alpha
    },

    #' @description
    #' Get the beta parameter
    #'
    #' @return Numeric. The background class (color) penalty/favor.
    beta = function(){
      private$.beta
    },

    #' @description
    #' Get the seeds
    #'
    #' @return A list of seed specifications.
    seeds = function(){
      private$.seeds
    }

  ),

  private = list(
    .alpha = NULL,
    .beta = NULL,
    .seeds = NULL,
    .ncolors = NULL
  )
)
