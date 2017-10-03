# TO DO: Fill out roxygen2 documentation.
#
#' @title SQP algorithm with partial QR for optimizing mixtures.
#' 
#' @description Sequential quadratic programming (SQP) method with
#'     partial QR-based approximate computation of gradients and Hessians
#'     for efficiently solving the Mixture Distribution Optimization
#'     Problem.
#' 
#' @export
#' @importFrom rjulia r2j
#' @importFrom rjulia j2r
#' @importFrom rjulia julia_void_eval
mixsqp <- function (L, x, convtol = 1e-8, pqrtol = 1e-8, eps = 1e-8,
                    sptol = 1e-3, maxiter = 100, maxqpiter = 100,
                    seed = 1, verbose = TRUE) {
    
  # Pass arguments from R to Julia.
  r2j(L,"L")
  r2j(x,"x")
  r2j(eps,"eps")
  r2j(sptol,"sptol")
  
  # Load the Julia code.
  julia_void_eval("include(\"../inst/code/mixsqp.jl\")")
  
  # Run the SQP algorithm.
  rjulia::julia_void_eval('out = mixsqp(L,x)');
  
  # Construct the return value.
  return(rjulia::j2r("out[1]"))
  # B <- rjulia::j2r('temp[2]'); # loglik at the solution
  # niter <- rjulia::j2r('temp[3]'); # number of 
  # converged <- rjulia::j2r('~temp[4]');
  # return(list(pihat = pihat, B=B, 
  #             niter = niter, converged=converged))
}

#' @importFrom stats rnorm
#' @importFrom stats rt
get_sample = function(n,seed=2017) {
  set.seed(seed)
  n1 = floor(n/2)
  z = c(rnorm(n1),4*rt(n-n1,df=6))
  z = z[order(abs(z))]
  return(z)
}

get_matrix_lik <- function(z,m = 1.1) {
  data = ashr::set_data(z,1)
  grid = ashr:::autoselect.mixsd(data, mult=m, mode=0, mixcompdist="normal", grange=c(-Inf,Inf))
  grid = c(0,grid)
  k = length(grid)
  g  = ashr::normalmix(rep(1/k,k),rep(0,k),grid)
  llik <- t(ashr:::log_comp_dens_conv.normalmix(g,data))
  L = llik - apply(llik, 1, max)
  L = exp(L)
  return(L)
}

mixopt_ex = function(n, m = 1.1) {
  cat("sample z from 0.5 * N(0,1) + 0.5 g\n");
  z = get_sample(10000);
  cat("compute likelihood : ");
  t1 = system.time(L <- get_matrix_lik(z,m = m));
  cat(t1[3],"sec\n");
  cat("convex programming : ");
  t2 = system.time(x <- mixopt(L));
  cat(t2[3],"sec\n");
  return(x)
}
