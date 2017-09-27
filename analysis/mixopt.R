## introduction to be written 

#'
mixopt_init = function(julia_package_install = FALSE){
  if (!require("rjulia",character.only = TRUE)){
    install.packages("rjulia",dep=TRUE)
    if(!require(rjulia,character.only = TRUE)) stop("Package not found")
  }
  cat("R package 'rjuila' installed\n")
  rjulia::julia_init();
  cat("Julia initiaized\n")
  rjulia::r2j(getwd(),"working_directory")
  rjulia::julia_void_eval('cd(working_directory[1])')
  cat("Working directory has changed in Julia\n")
  if(julia_package_install == TRUE){
    rjulia::julia_void_eval('Pkg.add("LowRankApprox")');
    cat("Julia package 'LowRankApprox' installed\n");
  }
}

#' @param : matrix_lik : likelihood matrix of dimension n times k
#'                 eps :  
#' @return pihat : 
#' 
#' 

get_sample = function(n,seed=2017){
  set.seed(seed)
  n1 = floor(n/2)
  z = c(rnorm(n1),4*rt(n-n1,df=6))
  z = z[order(abs(z))]
  return(z)
}

get_matrix_lik = function(z,m = 1.1){
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

mixopt = function(matrix_lik, eps = 1e-8, tol = 1e-8, sptol = 1e-3){
  
  # pass arguments from R to Julia
  rjulia::r2j(matrix_lik,"matrix_lik");
  rjulia::r2j(eps,"eps"); rjulia::r2j(tol,"tol"); rjulia::r2j(sptol,"sptol");
  
  # load Julia code
  rjulia::julia_void_eval('include("sqp.jl")');
  
  # run convex programming
  rjulia::julia_void_eval('temp = sqp(matrix_lik,eps[1],tol[1],sptol[1])');
  
  # get return value
  pihat <- rjulia::j2r('temp[1]'); # solution
  B <- rjulia::j2r('temp[2]'); # loglik at the solution
  niter <- rjulia::j2r('temp[3]'); # number of 
  converged <- rjulia::j2r('~temp[4]');
  return(list(pihat = pihat, B=B, 
              niter = niter, converged=converged))
}

mixopt_ex = function(n, m = 1.1){
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
