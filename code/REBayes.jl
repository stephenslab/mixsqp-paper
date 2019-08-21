# Fit mixture model using MOSEK interior-point solver, called by the
# KWDual function from the REBayes package.
function REBayes(L)
  @rput L;
  R"library(REBayes);
    n      <- nrow(L)
    k      <- ncol(L)
    timing <- system.time(out <- KWDual(L,rep(1,k),rep(1,n)/n));
    timing <- timing[3];
    x      <- out$f;
    x      <- x/sum(x);"
  @rget x;
  @rget timing;
  return x, timing
end
