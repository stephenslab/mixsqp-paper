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
List mixsqp_rcpp (const arma::mat& L, const arma::vec& x0,
		  double convtol, double pqrtol, double eps,
		  double sptol, int maxiter, int maxqpiter,
		  bool verbose) {

  // Get the number of rows (n) and columns (k) of the conditional
  // likelihood matrix.
  int n = L.n_rows;
  int k = L.n_cols;

  // Print a brief summary of the analysis, if requested.
  if (verbose) {
    Rprintf("Running SQP algorithm with the following settings:\n");
    Rprintf("- %d x %d data matrix\n",n,k);
    Rprintf("- convergence tolerance = %0.2e\n",convtol);
    Rprintf("- zero threshold        = %0.2e\n",sptol);
  }

  // PREPARE DATA STRUCTURES FOR OPTIMIZATION ALGORITHM
  // Initialize storage for the outputs obj, gmin, nnz and nqp.
  arma::vec obj(maxiter);
  arma::vec gmin(maxiter);
  arma::vec nnz(maxiter);
  arma::vec nqp(maxiter);
  
  // Initialize the solution.
  arma::vec x = x0;
  
  // Initialize storage for matrices and vectors used in the
  // computations below.
  arma::rowvec g(k);   // Vector of length k storing the gradient.
  arma::vec    u(n);   // Vector of length n storing L*x + eps or its log.
  arma::mat    H(k,k); // k x k matrix storing Hessian.
  arma::mat    Z(n,k); // n x k matrix storing  Z = diag(1/(L*x + eps))*L.
  arma::mat    I(k,k); // k x k diagonal matrix eps*I.
  arma::uvec   t(k);   // Temporary unsigned integer vector result of length k.

  // Initialize storage for additional matrices and vectors used for
  // the inner loop.
  uint       ns;         // Number of variables in subproblem.
  arma::uvec indmem(k);  // Storage for nonzero indices of x.
  arma::vec  ymem(k);    // Storage for solution to subproblem.
  arma::vec  gsmem(k);   // Storage for gradient of subproblem.
  arma::vec  dsmem(k);   // Storage for search direction in subproblem.
  arma::mat  Hsmem(k,k); // Storage for Hessian of subproblem.
  
  // This is used in computing the Hessian matrix.
  I  = arma::eye(k,k);
  I *= eps;
  
  // Initialize some loop variables used in the loops below.
  int j = 0;

  // Print the column labels for reporting the algorithm's progress.
  if (verbose)
    Rprintf("iter       objective -min(g+1) #nnz #qp\n");

  // Repeat until we reach the maximum number of outer loop iterations.
  for (int i = 0; i < maxiter; i++) {

    // COMPUTE GRADIENT AND HESSIAN
    // Compute Z = diag(1/(L*x + eps)) * L.
    u = L * x + eps;
    Z = L;
    Z.each_col() /= u;

    // Compute the gradient g = -Z'*1/n.
    g = arma::sum(Z,0);
    g = -g/n;

    // Compute the Hessian H = Z'*Z/n + eps*I.
    H = Z.t() * Z;
    H = H/n + I;

    // Report on the algorithm's progress. Here we compute: the value
    // of the objective at x (obj); the smallest gradient value
    // (gmin), which is used as a convergence criterion; the number of
    // nonzeros in the solution (nnz); and the number of inner-loop
    // iterations (nqp).
    //
    // TO DO: The L * x matrix operation here used to compute the
    // objective function could slow down the algorithm when number of
    // QR factors in the partial QR is much smaller than / k. We need
    // to think of a way to avoid this by having an option to not
    // output the objective function at each iteration, and/or make
    // sure that this objective function operation is not included in
    // the timing.
    //
    u       = log(u);
    t       = (x > sptol);
    obj[i]  = -sum(u);
    gmin[i] = 1 + min(g);
    ns      = sum(t);
    nnz[i]  = ns;
    nqp[i]  = j;
    if (verbose)
      Rprintf("%4d %0.8e %+0.2e %4d %3d\n",i,obj[i],-gmin[i],nnz[i],j);
      
    // Check convergence.
    // 
    // TO DO: Don't we want some sort of convergence tolerance
    // parameter here? e.g., min(g+1) >= d, where d is a positive
    // number near zero.
    //
    if (gmin[i] >= 0)
      break;

    // Initialize the solution to the QP subproblem (y).
    arma::vec y(ymem.memptr(),ns,false,true);
    y.fill(1/(double) ns);
      
    // Get the gradient and Hessian for the QP subproblem.
    arma::uvec ind(indmem.memptr(),k,false,true);
    arma::vec  gs(gsmem.memptr(),k,false,true);
    arma::vec  ds(gsmem.memptr(),k,false,true);
    arma::mat  Hs(Hsmem.memptr(),k,k,false,true);
    
    ind = find(x > sptol);
    gs  = g.elem(ind);
    Hs  = H.elem(ind,ind);
    
    // Run active set method to solve the QP subproblem.
    for (j = 0; j < maxqpiter; j++) {
          
      // Define the smaller QP subproblem.
      ds = Hs*y + 2*gs + 1;

      // Solve the smaller problem.
      // p      = sparse(zeros(k));
      // p_s    = -H_s\d_s;
      // p[ind] = p_s;

      // Check convergence.
      //  
      // TO DO: Why is convergence based on the norm of the search
      // direction? Please relate to KKT conditions.
    //   if (norm(p_s) < convtol) {
            
    // 	// Compute the Lagrange multiplier.
    //     lambda = d - minimum(d_s);
    //     if all(lambda .>= 0)
    //       break;
    //     else {
            
    //       // TO DO: Explain what ind and ind_min are for.
    //       ind_min = findmin(lambda)[2];
    //       ind     = sort([ind; ind_min]);
    //     }
    //   } else {
          
    //     // Find a feasible step length.
    //     alpha     = 1;
    //     alpha0    = -y[ind]./p_s;
    //     ind_block = find(p_s .< 0);
    //     alpha0    = alpha0[ind_block];
    //     if (~isempty(ind_block)) {
    //       v, t = findmin(alpha0);
    //       if (v < 1) {

    //         // Blocking constraint.
    //         ind_block = ind[ind_block[t]]; 
    //         alpha     = v;
              
    //         // Update working set if there is a blocking constraint.
    //         deleteat!(ind,find(ind - ind_block .== 0));
    //       }
    //     }
          
    // 	// Move to the new "inner loop" iterate (y) along the search
    //     // direction.
    //     y = y + alpha * p;
    //   }
    }
    
    // Update the solution.
    x = y;
  }
  
  return List::create(Named("g") = g,Named("H") = H);
}
