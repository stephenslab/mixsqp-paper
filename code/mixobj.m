

y <- drop(L %*% x) + e
 if (all(y > 0))
   return(-sum(w * log(y)))
