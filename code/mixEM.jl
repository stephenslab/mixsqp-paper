# Fit a mixture model using EM. Input argument L is the n x m
# conditional likelihood matrix, where n is the number of samples and
# m is the number of mixture components; optional input argument w is
# the initial estimate of the mixture weights.
function mixEM(L; w = ones(size(L,2))/size(L,2), maxiter = 10000,
               tol = 1e-6, eps = 1e-8)

  # Get the number of rows (n) and columns (k) of the conditional
  # likelihood matrix.
  n = size(L,1);
  k = size(L,2);

  # Initialize storage for outputs obj, maxd and timing.
  obj    = zeros(maxiter);
  maxd   = zeros(maxiter);
  timing = zeros(maxiter);
    
  # Compute the objective function value at the initial iterate.
  iter      = 1;
  obj[iter] = -sum(log.(L * w + eps));
    
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  for iter = 2:maxiter

    # Start timing the SQP iteration.
    tic();
      
    # Save the current estimate of the mixture weights.
    w0 = w;

    # E STEP
    # ------                       
    # Compute the posterior probabilities.
    P = L * spdiagm(w);
    P = P ./ repmat(sum(P,2) + eps,1,k);

    # M STEP
    # ------                       
    # Update the mixture weights.
    w = mean(P,1)'[:];
    
    # COMPUTE OBJECTIVE
    # -----------------
    obj[iter] = -sum(log.(L * w + eps));
      
    # CHECK CONVERGENCE
    # -----------------
    # Convergence is reached when the maximum difference between the
    # mixture weights at two successive iterations is less than the
    # specified tolerance, or when objective increases.
    maxd[iter] = maximum(abs.(w - w0));
    timing[iter] = toq();
    if maxd[iter] < tol
      break
    end
  end
      
  # Return the mixture weights and other optimization info.
  return w, obj[1:iter], maxd[1:iter], timing[1:iter]
end        

