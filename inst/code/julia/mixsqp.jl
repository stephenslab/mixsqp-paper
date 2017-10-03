using LowRankApprox

# TO DO: Briefly explain here what this function does, and how to use it.
# TO DO: Describe the function inputs and outputs.
# TO DO: Maybe find a better name for input argumnet "convtol"?
# TO DO: Add per-iteration timing info.
function mixsqp(L, x; convtol = 1e-8, pqrtol = 1e-8, eps = 1e-8,
                sptol = 1e-3, maxiter = 100, maxqpiter = 100,
                seed = 1, verbose = true)

  # Start timing the function.
  tic();
    
  # Get the number of rows (n) and columns (k) of the conditional
  # likelihood matrix.
  n = size(L,1);
  k = size(L,2);

  # Compute partial QR decomposition with relative precision "tol",
  # then retrieve the permutation matrix, P. For details on the
  # partial QR, see https://github.com/klho/LowRankApprox.jl.
  srand(1);
  F = pqrfact(L,rtol = pqrtol);
  P = sparse(F[:P]);

  # Summarize the analysis here.
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

  # Initialize storage for the outputs obj, ming, nnz and nqp.
  obj    = zeros(maxiter);
  ming   = zeros(maxiter);
  nnz    = zeros(maxiter);
  nqp    = zeros(maxiter);
  timing = zeros(maxiter);
    
  # Initialize loop variables used in the loops below so that they
  # are available outside the scope of the loop.
  i = 0;
  j = 0;
  D = 0;

  # Print the column labels for reporting the algorithm's progress.
  if verbose
    @printf("iter       objective -min(g+1) #nnz #qp\n")
  end

  # QP subproblem start.
  for i = 1:maxiter

    # Start timing the iteration.
    tic();
      
    # Compute the gradient and Hessian, using the partial QR
    # decomposition to increase the speed of these computations.
    #
    # TO DO: Maybe simplify this code by first setting Q = F[:Q] and 
    # R = F[:R]?
    #
    D = 1./(F[:Q]*(F[:R]*(P'*x)) + eps);
    g = -P * F[:R]' * (F[:Q]'*D)/n;
    H = P * F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q]) * F[:R] * P'/n +
        pqrtol * eye(k);

    # Report on the algorithm's progress.
    if verbose
      obj[i]  = -sum(log.(L * x + eps));
      ming[i] = minimum(g + 1);
      nnz[i]  = sum(x .> sptol);
      nqp[i]  = j;
      @printf("%4d %0.8e %+0.2e %4d %3d\n",i,obj[i],-ming[i],nnz[i],j);
    end
      
    # Check convergence.
    #  
    # TO DO: Don't we want some sort of convergence tolerance
    # parameter here? e.g., min(g+1) >= d, where d is a positive
    # number near zero.
    if minimum(g + 1) >= 0
      break
    end
      
    # Initialize the solution to the QP subproblem (y).
    ind    = find(x .> sptol);
    y      = sparse(zeros(k));
    y[ind] = 1/length(ind);

    # Run active set method to solve the QP subproblem.
    for j = 1:maxqpiter
          
      # Define the smaller QP subproblem.
      s   = length(ind);
      H_s = H[ind,ind];
      d   = H*y + 2*g + 1;
      d_s = d[ind];

      # Solve the smaller problem.
      p      = sparse(zeros(k));
      p_s    = -H_s\d_s;
      p[ind] = p_s;

      # Check convergence.
      #  
      # TO DO: Why is convergence based on the norm of the search
      # direction? Please relate to KKT conditions.
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
          v, t = findmin(alpha0);
          if v < 1

            # Blocking constraint.
            ind_block = ind[ind_block[t]]; 
            alpha     = v;
              
            # Update working set if there is a blocking constraint.
            deleteat!(ind,find(ind - ind_block .== 0));
          end
        end
          
        # Move to the new "inner loop" iterate (y) along the search
        # direction.
        y = y + alpha * p;
      end
    end

    # Update the solution to the original optimization problem.
    x = y;

    # Get the elapsed time for the ith iteration.
    timing[i] = toq();
  end

  # Return: (1) the solution (after zeroing out any values below the
  # tolerance); (2) the value of the objective at each iteration; (3)
  # the minimum gradient value of the modified objective at each
  # iteration; (4) the number of nonzero entries in the vector at each
  # iteration; and (5) the number of inner iterations taken to solve
  # the QP subproblem at each outer iteration.
  x[x .< sptol] = 0
  totaltime = toq();  
  if verbose
    @printf("Optimization took %d iterations and %0.4f seconds.\n",i,totaltime)
  end
  return Dict([("x",full(x)), ("totaltime",totaltime), ("obj",obj[1:i]),
               ("ming",ming[1:i]), ("nnz",nnz[1:i]), ("nqp",nqp[1:i]),
               ("timing",timing[1:i])])
end
