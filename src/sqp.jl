using LowRankApprox

# TO DO: Briefly explain here what this function does, and how to use it.
# TO DO: Maybe find a better name for input argumnet "convtol"?
function sqp(L, x; convtol = 1e-8, pqrtol = 1e-8, eps = 1e-8,
             sptol = 1e-3, maxiter = 100, maxqpiter = 100,
             verbose = true)

  # Get the number of rows (n) and columns (k) of the conditional
  # likelihood matrix.
  n = size(L,1);
  k = size(L,2);

  # Compute partial QR decomposition with relative precision "tol",
  # then retrieve the permutation matrix, P. For details on the
  # partial QR decomposition, see
  # https://github.com/klho/LowRankApprox.jl.
  F = pqrfact(L,rtol = pqrtol);
  P = sparse(F[:P]);

  # TO DO: Add summary of analysis here:
  if verbose
    p   = F[:p];
    err = maximum(abs.(F[:Q]*F[:R] - L[:,p]));
    @printf("Running SQP algorithm with the following settings:\n")
    @printf("- %d x %d data matrix\n",n,k)
    @printf("- convergence tolerance = %0.2e\n",convtol)
    @printf("- zero threshold        = %0.2e\n",sptol)  
    @printf("- partial QR tolerance  = %0.2e\n",pqrtol)
    @printf("- partial QR max. error = %0.2e\n",err)
  end

  # Initialize loop variables used in the loops below so that they
  # are available outside the scope of the loop.
  i = 0;
  D = 0;

  # Report the initial algorithm state.
  if verbose
    f = -sum(log.(L * x + eps));
    @printf("iter     objective\n")
    @printf("   0 %0.6e\n",f)
  end
    
  # QP subproblem start.
  for i = 1:maxiter
        
    # Compute the gradient and Hessian, using the partial QR
    # decomposition to increase the speed of these computations..
    D = 1./(F[:Q]*(F[:R]*(P'*x)) + eps);
    g = -P * F[:R]' * (F[:Q]'*D)/n;
    H = P * F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q]) * F[:R] * P'/n +
        pqrtol * eye(k);

    # Initialize solution to QP subproblem.
    ind    = find(x .> sptol);
    y      = sparse(zeros(k));
    y[ind] = 1/length(ind);

    # Run active set method to solve QP subproblem.
    for j = 1:maxqpiter
          
      # Define the smaller problem.
      s   = length(ind);
      H_s = H[ind,ind];
      d   = H*y + 2*g + 1;
      d_s = d[ind];

      # Solve the smaller problem.
      p      = sparse(zeros(k));
      p_s    = -H_s\d_s;
      p[ind] = p_s;

      # Check convergence.
      # TO DO: Why is convergence based on the norm of the search
      # direction? Relate to KKT conditions.
      if norm(p_s) < convtol
            
        # Compute the Lagrange multiplier.
        lambda = d - minimum(d_s);
        if all(lambda .>= 0)
          break;
        else
            
          # TO DO: Explain what ind and ind_min are here.
          ind_min = findmin(lambda)[2];
          ind     = sort([ind; ind_min]);
        end
      else
          
        # Find a feasible step length.
        alpha     = 1;
        alpha0    = -y[ind]./p_s;
        ind_block = find(p_s .< 0);
        alpha0    = alpha0[ind_block];
        if ~isempty(ind_block)
          t = findmin(alpha0);
          if t[1] < 1

            # Blocking constraint.
            ind_block = ind[ind_block[t[2]]]; 
            alpha     = t[1];
              
            # Update working set if there is a blocking constraint.
            deleteat!(ind,find(ind - ind_block .== 0));
          end
        end
          
        # Move to the new iterate (y) along the search direction.
        y = y + alpha * p;
      end
    end
    x = y;

    # Check convergence.
    if minimum(g + 1) >= 0
      break
    end
  end
  x[x .< sptol] = 0
  return full(x), sum(log.(D + eps)), i, i == maxiter
end
