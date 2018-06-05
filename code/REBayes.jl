function REBayes(L)
  @rput L;
  R"library(REBayes);
    n           <- nrow(L)
    k           <- ncol(L)
    timing      <- system.time(out <- KWDual(L,rep(1,k),rep(1,n)/n));
    timing      <- timing[3];
    x           <- out$f;
    x[x < 1e-3] <- 0;
    x           <- x/sum(x);"
  @rget x;
  @rget timing;
  return x, timing
end
