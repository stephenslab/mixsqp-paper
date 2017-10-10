using LowRankApprox

# TO DO: Briefly explain here in a few sentences what this function
#        does, and how to use it.
# TO DO: Describe the function inputs and outputs.
# TO DO: Maybe find a better name for input argumnet "convtol"?
function mixsqp(L; x = [], w = [], convtol = 1e-8, pqrtol = 0, eps = 1e-8,
                sptol = 1e-3, maxiter = 100, maxqpiter = 100,
                seed = 1, verbose = true)

  # Start timing the function.
  tic();

  # Get the number of rows (n) and columns (k) of the conditional
  # likelihood matrix.
  n = size(L,1);
  k = size(L,2);

  # If specified (i.e., if x ~= []), we initialize with that x
  # Otherwise let x = ones(k)/k
  if length(x) == 0
    x = ones(k)/k
  end

  # If specified (i.e., if w ~= []]), we keep that value of w but normalized
  # so that sum(w) = n
  # Otherwise let w = ones(n)
  if length(w) == 0
    w = ones(n)
  else
    w = w/sum(w) * n
  end

  # If requested (i.e., if pqrtol > 0), compute the partial QR
  # decomposition with relative precision "tol", then retrieve the
  # permutation matrix, P. For details on the partial QR, see
  # https://github.com/klho/LowRankApprox.jl.
  if pqrtol > 0
    srand(1);
    F = pqrfact(L,rtol = pqrtol);
    P = sparse(F[:P]);
    Q = F[:Q];
    R = F[:R]
  end

  # Summarize the analysis here.
  if verbose
    @printf("Running SQP algorithm with the following settings:\n")
    @printf("- %d x %d data matrix\n",n,k)
    @printf("- convergence tolerance = %0.2e\n",convtol)
    @printf("- zero threshold        = %0.2e\n",sptol)
    if pqrtol > 0
      p   = F[:p];
      err = maximum(abs.(F[:Q]*F[:R] - L[:,p]));
      @printf("- partial QR tolerance  = %0.2e\n",pqrtol)
      @printf("- partial QR max. error = %0.2e\n",err)
    else
      @printf("- Exact derivative computation (partial QR not used).\n")
    end
  end

  # Initialize storage for the outputs obj, gmin, nnz and nqp.
  obj    = zeros(maxiter);
  gmin   = zeros(maxiter);
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

    # Compute the gradient and Hessian, optionally using the partial
    # QR decomposition to increase the speed of these computations.
    if pqrtol > 0
      D = sqrt.(w)./(Q*(R*(P'*x)) + eps);
      g = -P * R' * (Q'*(D.* sqrt.(w)))/n;
      H = P * R' * (Q'*Diagonal(D.^2)*Q) * R * P'/n + pqrtol * eye(k);
    else
      D = sqrt.(w)./(L*x + eps);
      g = -(L'*(D.*sqrt.(w)))/n;
      H = L' * Diagonal(D.^2) * L/n + eps * speye(k);
    end

    # Report on the algorithm's progress.
    #
    # TO DO: The L * x matrix operation here used to compute the
    # objective function could dramatically slow down the algorithm
    # when number of QR factors in the partial QR is much smaller than
    # k. We need to think of a way to avoid this by having an option
    # to not output the objective function at each iteration, and/or
    # make sure that this objective function operation is not included
    # in the timing.
    #
    obj[i]  = -sum(w .* log.(L * x + eps));
    gmin[i] = minimum(g + 1);
    nnz[i]  = sum(x .> sptol);
    nqp[i]  = j;
    if verbose
      @printf("%4d %0.8e %+0.2e %4d %3d\n",i,obj[i],-gmin[i],nnz[i],j);
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

          # TO DO: Explain what ind and ind_min are for.
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
  x = x/sum(x);
  totaltime = toq();
  if verbose
    @printf("Optimization took %d iterations and %0.4f seconds.\n",i,totaltime)
  end
  return Dict([("x",full(x)), ("totaltime",totaltime), ("obj",obj[1:i]),
               ("gmin",gmin[1:i]), ("nnz",nnz[1:i]), ("nqp",nqp[1:i]),
               ("timing",timing[1:i])])
end
