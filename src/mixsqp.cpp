#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

using namespace Rcpp;

// SQP algorithm for optimizing mixtures. For more information, see
// the help and comments accompanying the "mixsqp" function in R.
// 
// [[Rcpp::export]]
arma::vec mixsqp_rcpp (const arma::mat& L, const arma::vec& x0,
		       double convtol, double pqrtol, double eps, 
		       double sptol, int maxiter, int maxqpiter, 
		       bool verbose) {

  // Get the number of rows (n) and columns (k) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int k = L.n_cols;

  Rprintf("- %d x %d data matrix\n",n,k);

// - convergence tolerance = 1.00e-08
// - zero threshold        = 1.00e-03
// - partial QR tolerance  = 1.00e-08
// - partial QR max. error = 4.63e-08

  arma::vec x = x0;
  return x;
}
