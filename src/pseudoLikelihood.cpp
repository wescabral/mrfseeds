#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

 //' Compute Log Pseudo-Likelihood
 //'
 //' Calculates the log pseudo-likelihood for the MRF model given observations
 //' and parameters.
 //'
 //' @param z Integer matrix of hidden field values
 //' @param alpha Double - uniformity parameter
 //' @param beta Double - global background color penalty/favor
 //' @param ncolors Integer - number of colors
 //' @param distMatrix Numeric matrix - pre-calculated distances (n_pixels x n_seeds)
 //' @param seedColors Integer vector - color associated with each seed
 //' @param seedDeltas Numeric vector - delta parameter for radial attraction/repulsion of each seed
 //'
 //' @return Double - log pseudo-likelihood value
 //'
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]

double computeLogPL(IntegerMatrix z,
                    double alpha,
                    double beta,
                    int ncolors,
                    NumericMatrix distMatrix,  // n_pixels x n_seeds, pré-calculada
                    IntegerVector seedColors,
                    NumericVector seedDeltas) {
  int nrows = z.nrow();
  int ncols = z.ncol();
  int n = nrows * ncols;
  int nseeds = seedColors.size();

  double logPL = 0.0;

  for (int pos = 0; pos < n; pos++) {
    int pos_x = pos % nrows;
    int pos_y = pos / nrows;

    int observed = z(pos_x, pos_y);

    std::vector<int> neighbors;
    neighbors.reserve(4);
    if (pos_x > 0)         neighbors.push_back(z(pos_x - 1, pos_y));
    if (pos_x < nrows - 1) neighbors.push_back(z(pos_x + 1, pos_y));
    if (pos_y > 0)         neighbors.push_back(z(pos_x, pos_y - 1));
    if (pos_y < ncols - 1) neighbors.push_back(z(pos_x, pos_y + 1));

    double max_energy = R_NegInf;
    int observed_idx = 0;
    NumericVector energies(ncolors);

    for (int v = 0; v < ncolors; v++) {
      double contrib = 0.0;
      for (size_t j = 0; j < neighbors.size(); j++) {
        contrib += (v == neighbors[j]) ? alpha : 0.0;
      }

      double energy = contrib + (v == 0 ? beta : 0.0);

      for (int s = 0; s < nseeds; s++) {
        if (v == seedColors[s]) {
          energy += seedDeltas[s] / (1.0 + distMatrix(pos, s));
        }
      }

      energies[v] = energy;
      if (energy > max_energy) max_energy = energy;
      if (v == observed) observed_idx = v;
    }

    // truque de subtrair o máximo
    double sum_exp = 0.0;
    for (int v = 0; v < ncolors; v++) {
      sum_exp += std::exp(energies[v] - max_energy);
    }

    logPL += energies[observed_idx] - max_energy - std::log(sum_exp);
  }

  return logPL;
}
