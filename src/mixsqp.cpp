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
Rcpp::List mixsqp_rcpp (const arma::mat& L, const arma::vec& x0,
			double convtol, double pqrtol, double eps,
			double sptol, int maxiter, int maxqpiter,
			bool verbose) {

  // Get the number of rows (n) and columns (k) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int k = L.n_cols;

  // Print a brief summary of the analysis.
  Rprintf("Running SQP algorithm with the following settings:\n");
  Rprintf("- %d x %d data matrix\n",n,k);
  Rprintf("- convergence tolerance = %0.2e\n",convtol);
  Rprintf("- zero threshold        = %0.2e\n",sptol);

  // Initialize the solution.
  arma::vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // optimization algorithm below.
  arma::rowvec g(k);
  arma::vec u(n);
  arma::mat U(n,n);
  arma::mat H(k,k);
  arma::mat Z(n,k);
  arma::mat I(k,k);
  I = eps * arma::eye(k,k);
  
  // Repeat until we reach the maximum number of outer loop iterations.
  for (int i = 0; i < maxiter; i++) {

    // Compute the gradient.
    u = L * x + eps;
    u = 1 / u;
    U = arma::diagmat(u);
    Z = U * L;
    g = arma::sum(Z,0);
    g = -g / n;

    // Compute the Hessian.
    H = Z.t() * Z;
    H = H/n;
    H = H + I;
  }
  
  return List::create(Named("g") = g,Named("H") = H);
}
