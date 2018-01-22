#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
// 
// [[Rcpp::depends(RcppArmadillo)]]

// Some functions and classes uses frequently in the code below.
using arma::vec;

arma::vec mixsqp_rcpp
