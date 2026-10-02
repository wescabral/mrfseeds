# mrfseeds

Image generation and bayesian parameter inference package for MRFs with centrality effects

## Features
- Generate synthetic images using MRF models
- Bayesian parameter inference with Metropolis-within-Gibbs
- Seeds visualization

## Installation
```r
# Install from GitHub
devtools::install_github("wescabral/mrfseeds")
```

## Example
```r
library(mrfseeds)

# Generate a random image
real_seeds <- list(
  list(position = c(50, 50),   color = 1, delta = 20),
  list(position = c(100, 100), color = 1, delta = 20)
)
mod <- seedModel$new(alpha = 0, beta = 1, seeds = real_seeds, ncolors = 2)
img <- mod$sampleImage(dim = c(150, 150), steps = 80)
img$plot()

# Noisy image
Y <- img$matrix + rnorm(prod(dim(img$matrix)), sd = 0.1)
img$setMatrix(Y, continuous = TRUE)
img$plot()

# Infer parameters
priors <- list(
  shape_alpha = 2, rate_alpha = 1, 
  shape_beta  = 10, rate_beta  = 0.5,
  shape_delta = 5, rate_delta = 0.5,
  m = 0, tau = 0.001,
  a = 10, b = 10
)
sampler <- seedBayesianCont$new(image = img, seeds_info = seeds_info, ncolors = 2, priors = priors)
sampler$plotZ()   # First look at the cleansed image

# init = NULL would initialize over PL inference
res <- sampler$run(n_iter = 10000, step_size = 0.01, init = list(alpha = 1, beta = 2, deltas = c(20, 20, 20)))


# Generate results
res$plot_logpl_trace()          # Convergence analysis
res$plot_image(n_samples = 300, burn_in = 3000) # Circles for influence area estimates
res$plot_traces(burn_in = 3000)
res$plot_densities(burn_in = 3000)
print(res, burn_in = 3000)
res$posterior_mean(burn_in = 3000)
res$credible_interval(level = 0.95, burn_in = 3000)
res$plot_positions()
```