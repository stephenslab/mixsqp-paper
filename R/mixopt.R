#' @title SQP algorithm with partial QR for optimizing mixtures.
#' 
#' @description Sequential quadratic programming (SQP) method with
#'     partial QR-based approximate computation of gradients and Hessians
#'     for efficiently solving the mixture distribution optimization
#'     problem.
#'
#' @param L \eqn{n x k} data matrix with positive entries. For fitting
#'     mixture models, n is the number of samples and k is the number of
#'     mixture components, and the entries of L are the conditional
#'     likelihoods for each sample (row) and each mixture component
#'     (column).
#'
#' @param x Vector of length k containing the initial estimate of the
#'     solution. For fitting mixture models, x is the initial estimate of
#'     the mixture weights.
#'
#' @param convtol Tolerance for stopping criterion in SQP method.
#'    TO DO: explain more precisely what the stopping criterion is.
#'
#' @param pqrtol Relative precision used to determine number of
#'    factors in QR decomposition of matrix L. See the LowRankApprox
#'    Julia package for details.
#'
#' @param eps A small positive number used to stabilize some of the
#'    numerical operations.
#'
#' @param sptol All entries of x below this value are treated as zero.
#'
#' @param maxiter Maximum number of outer loop iterations.
#'
#' @param maxqpiter Maximum number of inner loop iterations for
#'     solving quadratic subproblem.
#'
#' @param seed Seed for pseudorandom number generator passed to
#'     function \code{srand} in Julia.
#'
#' @param verbose If \code{verbose = TRUE}, print progress of algorithm
#'     to console.
#'
#' @return \code{mixsqp} returns a list containing the following components:
#'     \item{x}{Estimated solution.}
#'     \item{totaltime}{Total elapsed time to compute solution.}
#'     \item{obj}{Value of objective function at each iteration of the
#'                SQP algorithm.}
#'     \item{gmin}{Minimum value of the modified objective gradient,
#'                 which is used to assess convergence.}
#'     \item{nnz}{Number of nonzeros (below \code{sptol}) in solution at
#'                each iteration.}
#'     \item{timing}{Elapsed time for each SQP iteration.}
#' 
#' @details
#'     TO DO: Details about the algorithm go here, such as a
#'     mathematical description of the optimization problem, and a brief
#'     description of how the SQP algorithm works.
#'
#' @examples
#' # Fit mixture model to "normmix" data set, using partial QR
#' # decomposition to speed up computation.
#' data(normmix.data)
#' L   <- normmix.data$L
#' out <- mixsqp(L,pqrtol = 1e-8)
#' cat("Compare SQP solution against the IP solution:\n")
#' print(round(data.frame(ip = normmix.data$w,sqp = out$x),digits = 6))
#' 
#' @export
#' @importFrom rjulia r2j
#' @importFrom rjulia j2r
#' @importFrom rjulia jDo
mixsqp <- function (L, x, convtol = 1e-8, pqrtol = 0, eps = 1e-8,
                    sptol = 1e-3, maxiter = 100, maxqpiter = 100,
                    seed = 1, verbose = TRUE) {

  # Get the number of mixture components.
  k <- ncol(L)  
    
  # Set default initial estimate for x.
  if (missing(x))
    x <- rep(1/k,k)
    
  # Pass arguments from R to Julia. Since r2j treats everything as a
  # matrix, additional (somewhat painful) steps need to be taken to
  # pass scalars, including integers and logicals.
  r2j(L,"L")
  r2j(x,"x")
  r2j(convtol,"convtol")
  r2j(pqrtol,"pqrtol")
  r2j(eps,"eps")
  r2j(sptol,"sptol")
  r2j(maxiter,"maxiter")
  r2j(maxqpiter,"maxqpiter")
  r2j(seed,"seed")
  r2j(verbose,"verbose")
  jDo("convtol   = convtol[1];")
  jDo("pqrtol    = pqrtol[1];")
  jDo("eps       = eps[1];")
  jDo("sptol     = sptol[1];")
  jDo("maxiter   = Integer(maxiter[1]);")
  jDo("maxqpiter = Integer(maxqpiter[1]);")
  jDo("seed      = Integer(seed[1]);")
  jDo("verbose   = verbose[1];")
  
  # Load the Julia code. Note that this code may not work on all
  # systems (e.g., Windows)---in the future it would be better to
  # create a separate mixopt Julia package so that "include" calls
  # inside the "jDo" call  are not needed.
  mixsqp.file <- system.file("code/julia/mixsqp.jl",package = "mixopt")
  jDo(sprintf("include(\"%s\")",mixsqp.file))
  
  # Run the SQP algorithm.
  jDo(paste("out = mixsqp(L,x,convtol = convtol,pqrtol = pqrtol",
            "eps = eps,sptol = sptol,maxiter = maxiter",
            "maxqpiter = maxqpiter,seed = seed,verbose = verbose)",
            sep = ","))

  # Construct the return value. It seems (though not verified) that
  # assigning the individual "dictionary" (i.e., list) entries to
  # variables within Julia helps to prevent errors about memory not
  # being correctly mapped.
  #
  # TO DO: For convenience (e.g., to make it easier to draw plots),
  # the per-iteration information should be returned in a data frame.
  # 
  return(list(x         = j2r("x         = out[\"x\"];"),
              totaltime = j2r("totaltime = out[\"totaltime\"];"),
              obj       = j2r("obj       = out[\"obj\"];"),
              gmin      = j2r("gmin      = out[\"gmin\"];"),
              nnz       = j2r("nqp       = out[\"nqp\"];"),
              timing    = j2r("timing    = out[\"timing\"];")))
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
  t2 = system.time(x <- mixsqp(L));
  cat(t2[3],"sec\n");
  return(x)
}
