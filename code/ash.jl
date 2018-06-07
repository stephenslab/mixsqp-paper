# Basic Julia implementation of adaptive shrinkage with normal
# likelihood and mixture-of-centered-normals prior. The outputs are
# the fitted mixture weights (x), the standard deviations of the
# normal mixture components (sd), the posterior mean effect estimates
# (post_mean), and the runtime of the likelihood computations
# (t_likelihood), model fitting (t_fit) and posterior computations
# (t_posterior).
function ash(x, s; method = "mixSQP", gridmult = 1.3, lowrank = "svd",
             pqrtol = 1e-8)

  # Select the variances of the normal mixture components, then
  # compute the n x m conditional likelihood matrix.
  out, t_likelihood, bytes, gctime, memallocs = @timed if true
    sd = autoselectmixsd(x,s,gridmult = gridmult);
    L  = normlikmatrix(x,s,sd = sd);
  end
    
  # Fit the model using either the SQP or interior point solvers.
  if method == "mixSQP"
    fit, t_fit, bytes, gctime, memallocs = @timed mixSQP(L,pqrtol = pqrtol,
      lowrank = lowrank,eps = 1e-6,verbose = false);
    w = fit["x"];
  elseif method == "REBayes"
    w, t_fit = REBayes(L);
  end
    
  # Posterior calculations.
  out, t_posterior, bytes, gctime, memallocs = @timed if true
    ind = find(w .> 0);
    ps2 = sd[ind].^2;
    s2  = s.^2;
    t   = s2 .+ ps2';
    comp_post_mean  = (x * ps2') ./ t;
    comp_post_sd2   = (s2 * ps2') ./ t;
    comp_post_mean2 = comp_post_sd2 + comp_post_mean.^2;
    comp_post_prob  = L[:,ind] .* w[ind]';
    comp_post_prob  = comp_post_prob ./ sum(comp_post_prob,2);
    post_mean       = sum(comp_post_prob .* comp_post_mean,2);
    post_mean2      = sum(comp_post_prob .* comp_post_mean2,2);
  end
    
  return Dict([("x",x), ("timing-likelihood",t_likelihood),
               ("timing-fit",t_fit), ("timing-posterior",t_posterior)])
end
