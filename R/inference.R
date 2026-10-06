#' Seed Fit Result Class
#'
#' An R6 class that stores and visualizes the results of fitting a seed model
#' to image data using pseudo-likelihood maximization.
#'
#' @details
#' The `seedFitResult` class encapsulates:
#' - The fitted model with estimated parameters
#' - The log pseudo-likelihood value at convergence
#' - Convergence status
#' - Parameter estimates (alpha, beta, deltas)
#' - The original image used for fitting
#'
#' @examples
#' \dontrun{
#' # Fit a model (see seedModelFitter)
#' fitter <- seedModelFitter$new(image = img, seeds_info = seeds_info, ncolors = 2)
#' result <- fitter$fit()
#'
#' # Print summary
#' print(result)
#'
#' # Plot results with seed influence circles
#' result$plot()
#' }
#'
#' @export
seedFitResult <- R6::R6Class(
  "seedFitResult",
  public = list(

    #' @description
    #' Initialize a seed fit result
    #'
    #' @param model A `seedModel` object with fitted parameters
    #' @param logPL Numeric. The log pseudo-likelihood value at convergence
    #' @param convergence Integer. Convergence code from optimization routine
    #'   (0 = successful convergence, non-zero = failure)
    #' @param estimates A named numeric vector with parameter estimates:
    #'   alpha, beta, delta1, delta2, ...
    #' @param image An `imageData` object
    #'
    #' @return A seedFitResult object
    initialize = function(model, logPL, convergence, estimates, image) {
      private$.model       <- model
      private$.logPL       <- logPL
      private$.convergence <- convergence
      private$.estimates   <- estimates
      private$.image       <- image
      private$.seeds       <- model$seeds
    },

    #' @description
    #' Print a summary of the fit results
    #'
    #' @details
    #' Displays the convergence status, log pseudo-likelihood value,
    #' and estimated parameters
    print = function(...) {
      status <- if (private$.convergence == 0) "converged" else paste0("failed (code ", private$.convergence, ")")
      cat(sprintf("<seedFitResult> [%s]\n", status))
      cat(sprintf("  log-PL    : %.4f\n", private$.logPL))
      cat("  estimates :\n")
      for (nm in names(private$.estimates)) {
        val <- private$.estimates[[nm]]
        if (is.list(val)) {
          cat(sprintf("    %-10s  [%d seeds]\n", nm, length(val)))
        } else {
          cat(sprintf("    %-10s  %.4f\n", nm, val))
        }
      }
      invisible(self)
    },

    #' @description
    #' Plot the image with seed locations and influence area
    #'
    #' @details
    #' Overlays:
    #' - The original image
    #' - White crosses (+) marking seed positions
    #' - White circles showing the influence radius of each seed
    #'   (radius = delta/beta - 1, only if positive)
    #'
    #' @return A ggplot2 object
    plot = function(...) {
      beta  <- private$.model$beta
      seeds <- private$.model$seeds
      theta <- seq(0, 2 * pi, length.out = 300)

      p <- private$.image$plot() + ggplot2::coord_fixed()

      for (i in seq_along(seeds)) {
        seed  <- seeds[[i]]
        px    <- seed$position[1]
        py    <- seed$position[2]
        delta <- seed$delta

        p <- p + ggplot2::annotate("point", x = px, y = py,
                          shape = 3, size = 3, color = "white", stroke = 1)

        r <- if (beta > 0) delta / beta - 1 else NA_real_
        if (!is.na(r) && r > 0) {
          p <- p + ggplot2::annotate("path",
                            x = px + r * cos(theta),
                            y = py + r * sin(theta),
                            color = "white", linewidth = 0.7)
        }
      }

      p
    }
  ),

  active = list(

    #' @description
    #' Get the fitted model
    #'
    #' @return The `seedModel` object with fitted parameters
    model = function() private$.model,

    #' @description
    #' Get the log pseudo-likelihood at convergence
    #'
    #' @return Numeric. The log-PL value.
    logPL = function() private$.logPL,

    #' @description
    #' Get the convergence status
    #'
    #' @return Integer. 0 = converged, non-zero = failed
    convergence = function() private$.convergence,

    #' @description
    #' Get the parameter estimates
    #'
    #' @return A named numeric vector with alpha, beta, delta1, delta2, ...
    estimates = function() private$.estimates,

    #' @description
    #' Get the image used for fitting
    #'
    #' @return The `imageData` object
    image = function() private$.image,

    #' @description
    #' Get the fitted seeds
    #'
    #' @return A list of seed specifications with fitted delta values
    seeds = function() private$.seeds
  ),

  private = list(
    .model       = NULL,
    .logPL       = NULL,
    .convergence = NULL,
    .estimates   = NULL,
    .image       = NULL,
    .seeds       = NULL
  )
)

#' Seed Model Fitter Class
#'
#' An R6 class for fitting a seed model to image data via pseudo-likelihood maximization.
#'
#' @details
#' The `seedModelFitter` class:
#' - Initializes seed positions using k-means clustering
#' - Computes distance matrices for efficient likelihood evaluation
#' - Performs optimization to maximize the pseudo-likelihood
#' - Uses log-reparameterization for unconstrained optimization of positive parameters
#'
#' @examples
#' \dontrun{
#' # Create fitter
#' fitter <- seedModelFitter$new(
#'   image = img,
#'   seeds_info = list(list(color = 1, n_seeds = 2)),
#'   ncolors = 2
#' )
#'
#' # Fit the model
#' result <- fitter$fit(alpha_init = 1.0, beta_init = 0.5)
#' print(result)
#' }
#'
#' @export
seedModelFitter <- R6::R6Class(
  "seedModelFitter",
  public = list(

    #' @description
    #' Initialize a seed model fitter
    #'
    #' @param image An `imageData` object or a numeric matrix
    #' @param seeds_info A list of seed configuration objects. Each object should have:
    #'   - `color`: integer, the color/category to search for seeds
    #'   - `n_seeds`: integer, number of seeds to initialize for this color
    #' @param ncolors Integer. The number of color categories in the image
    #'
    #' @return A seedModelFitter object
    initialize = function(image, seeds_info, ncolors) {
      if(is.matrix(image)){
        image <- imageData$new(image)
      }
      private$.image    <- image
      private$.seeds    <- self$seeds_kmeans(seeds_info)
      private$.ncolors  <- ncolors
      nseeds     <- length(private$.seeds)
      seedMatrix <- do.call(rbind, lapply(private$.seeds, function(s) s$position))
      private$.seedMatrix <- seedMatrix
      private$.seedColors <- as.integer(sapply(private$.seeds, function(s) s$color))

      # Distance matrix precalc: n_pixels x n_seeds
      nrows <- image$dim[1]
      ncols <- image$dim[2]
      n     <- nrows * ncols

      pos_x <- (seq_len(n) - 1L) %% nrows
      pos_y <- (seq_len(n) - 1L) %/% nrows

      distMatrix <- matrix(0.0, nrow = n, ncol = nseeds)
      for (s in seq_len(nseeds)) {
        seed_x <- seedMatrix[s, 1] - 1
        seed_y <- seedMatrix[s, 2] - 1
        distMatrix[, s] <- sqrt((pos_x - seed_x)^2 + (pos_y - seed_y)^2)
      }
      private$.zMatrix    <- matrix(as.integer(image$matrix), nrow = nrows, ncol = ncols)
      private$.distMatrix <- distMatrix
    },

    #' @description
    #' Compute log pseudo-likelihood for given parameters
    #'
    #' @param alpha Numeric. The uniformity parameter
    #' @param beta Numeric. The background class (color) penalty/favor
    #' @param deltas Numeric vector. The influence radii for each seed (must be > 0)
    #'
    #' @return The log pseudo-likelihood value
    logPL = function(alpha, beta, deltas) {
      computeLogPL(
        z          = private$.zMatrix,
        alpha      = alpha,
        beta       = beta,
        ncolors    = private$.ncolors,
        distMatrix = private$.distMatrix,
        seedColors = private$.seedColors,
        seedDeltas = as.numeric(deltas)
      )
    },

    #' @description
    #' Initialize seed positions using k-means clustering
    #'
    #' @param seeds_info A list of seed configuration objects with `color` and `n_seeds`
    #'
    #' @return A list of initialized seed positions
    seeds_kmeans = function(seeds_info) {
      all_seeds <- list()

      for (config in seeds_info) {
        color <- config$color
        n_seeds <- config$n_seeds

        color_coords <- which(private$.image$matrix == color, arr.ind = TRUE)

        if (nrow(color_coords) < n_seeds) {
          stop(sprintf("Not enough color %d pixels (%d) for %d seeds. Need at least %d pixels.",
                       color, nrow(color_coords), n_seeds, n_seeds))
        }

        # Média pro caso de uma seed, kmeans pra mais de uma
        if (n_seeds == 1) {
          centroid <- colMeans(color_coords)
          centroids <- matrix(round(centroid, 0), nrow = 1)
        } else {
          kmeans_result <- kmeans(color_coords, centers = n_seeds, iter.max = 100)
          centroids <- round(kmeans_result$centers, 0)

          # Ordena as seeds
          distances_from_origin <- rowSums(centroids^2)
          centroid_order <- order(distances_from_origin)
          centroids <- centroids[centroid_order, ]
        }

        for (i in seq_len(n_seeds)) {
          all_seeds <- c(all_seeds, list(list(
            position = as.integer(centroids[i, ]),
            color = as.integer(color)
          )))
        }
      }

      return(all_seeds)
    },

    #' @description
    #' Fit the model by maximizing the pseudo-likelihood
    #'
    #' @details
    #' Uses the BFGS optimization algorithm with log-reparameterization for
    #' positive parameters (alpha and deltas) to ensure valid parameter values.
    #'
    #' @param alpha_init Numeric. Initial value for alpha. Default: 1.0
    #' @param beta_init Numeric. Initial value for beta. Default: 0.0
    #' @param deltas_init Numeric vector. Initial values for deltas.
    #'   If `NULL`, all initialized to 1.0
    #'
    #' @return A `seedFitResult` object with fitted parameters and diagnostics
    fit = function(alpha_init = 1.0, beta_init = 0.0, deltas_init = NULL) {
      nseeds <- length(private$.seeds)
      if (is.null(deltas_init)) deltas_init <- rep(1.0, nseeds)

      par_init <- c(log(alpha_init), beta_init, log(deltas_init))

      neg_logPL <- function(par) {
        alpha  <- exp(par[1])
        beta   <- par[2]
        deltas <- exp(par[seq(3, 2 + nseeds)])
        -self$logPL(alpha, beta, deltas)
      }

      result <- optim(par_init, neg_logPL, method = "BFGS")

      alpha_hat  <- exp(result$par[1])
      beta_hat   <- result$par[2]
      deltas_hat <- exp(result$par[seq(3, 2 + nseeds)])

      fitted_seeds <- mapply(
        function(s, d) list(position = s$position, color = s$color, delta = d),
        private$.seeds, deltas_hat,
        SIMPLIFY = FALSE
      )

      fitted_model <- seedModel$new(
        alpha   = alpha_hat,
        beta    = beta_hat,
        seeds   = fitted_seeds,
        ncolors = private$.ncolors
      )

      estimates <- c(
        alpha = alpha_hat,
        beta  = beta_hat,
        setNames(deltas_hat, paste0("delta", seq_along(deltas_hat)))
      )

      seedFitResult$new(
        model       = fitted_model,
        logPL       = -result$value,
        convergence = result$convergence,
        estimates   = estimates,
        image       = private$.image
      )
    },

    #' @description
    #' Print a summary of the fitter configuration
    print = function(...) {
      nseeds <- length(private$.seeds)
      dim    <- private$.image$dim
      cat("<seedModelFitter>\n")
      cat(sprintf("  image   : %d x %d  (%d colors)\n", dim[1], dim[2], private$.ncolors))
      cat(sprintf("  seeds   : %d\n", nseeds))
      for (i in seq_len(nseeds)) {
        s <- private$.seeds[[i]]
        cat(sprintf("    [%d] position = (%d, %d),  color = %d\n",
                    i, s$position[1], s$position[2], s$color))
      }
      invisible(self)
    }
  ),

  private = list(
    .image      = NULL,
    .seeds      = NULL,
    .ncolors    = NULL,
    .zMatrix    = NULL,
    .seedMatrix = NULL,
    .seedColors = NULL,
    .distMatrix = NULL
  )
)
