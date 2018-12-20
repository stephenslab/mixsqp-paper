# Input argument x specifies the initial iterate of the optimization
# algorithm. When x = -1, or when any of the entries of x are
# negative, the default setting for x is used.
#
# Since the MOSEK solver for the quadratic subproblem sometimes does
# not tolerate iterates in which most of the entries are nonzero
# ("dense" vectors), the default initial estimate for qpsubprob =
# "mosek" is a sparse vector in which only the first two entries are
# nonzero.
function mixsqp(L; x = -1,
                convergence_tol = 1e-8, eps_tol = 1e-8, sparse_tol = 1e-6,
                lowrank_method = "svd", lowrank_criterion = "by_rank",
                lowrank_rtol = 1e-8, lowrank_rank = minimum(size(L)),
                maxiter = 100, maxqpiter = 100,
                nullprior = 0, gradient_calc = "approx",
                linesearch = true, verbose = true)
    
# Get the number of rows (n) and columns (k) of L.
  n = size(L,1);
  k = size(L,2);

  # If any entries of input x are negative, set the initial estimate
  # to the default setting.
  #
  # When the MOSEK solver is used, the initial estimate is set to a
  # sparse vector in which all but two of the entries are zero.
  if any(x .< 0)
    x = ones(k)/k;
  end
    
  # If requested (i.e., if lowrank_rtol > 0), compute the partial QR
  # decomposition with relative precision "tol", then retrieve the
  # permutation matrix, P. For details on the partial QR, see
  # https://github.com/klho/LowRankApprox.jl.
    
  # Start timing for low-rank approximation of L.

  tic();
    
    if lowrank_method == "qr"
      if lowrank_criterion == "by_rank"
        F = pqrfact(L, rank = lowrank_rank);
        P = sparse(F[:P]);
      elseif lowrank_criterion == "by_rtol"
        F = pqrfact(L, rtol = lowrank_rtol);
        P = sparse(F[:P]);
      else
        error("lowrank_criterion must be either by_rank or by_rtol");
      end
          
    elseif lowrank_method == "svd"
      if lowrank_criterion == "by_rank"
        F = psvdfact(L, rank = lowrank_rank);
        S = Diagonal(F[:S]);
      elseif lowrank_criterion == "by_rtol"
        F = psvdfact(L, rtol = lowrank_rtol);
        S = Diagonal(F[:S]);
      else
        error("lowrank_criterion must be either by_rank or by_rtol");
      end
    end

  lowranktime = toq();
    
  # Summarize the analysis here.
  if verbose
    @printf("Running SQP algorithm with the following settings:\n")
    @printf("- %d x %d data matrix\n",n,k)
    @printf("- convergence tolerance = %0.2e\n",convergence_tol)
    @printf("- zero threshold        = %0.2e\n",sparse_tol)
    if lowrank_method == "qr"
      err = maximum(abs.(F[:Q]*F[:R]*P' - L));
      @printf("- partial QR tolerance  = %0.2e\n",lowrank_rtol)
      @printf("- partial QR max. error = %0.2e\n",err)
    elseif lowrank_method == "svd"
      err = maximum(abs.(F[:U]*S*F[:Vt] - L));
      @printf("- partial SVD tolerance  = %0.2e\n",lowrank_rtol)
      @printf("- partial SVD max. error = %0.2e\n",err)
    else
      @printf("- Exact derivative computation (partial QR not used).\n")
    end
  end

  # Initialize storage for the outputs obj, gmin, nnz and nqp.
  obj      = zeros(maxiter);
  gmin     = zeros(maxiter);
  nnz      = zeros(Int,maxiter);
  nqp      = zeros(Int,maxiter);
  timing   = zeros(maxiter);
  qptiming = zeros(maxiter);
    
  # Initialize loop variables used in the loops below so that they
  # are available outside the scope of the loop.
  i     = 0;
  j     = 0;
  D     = 0;
  t     = 0;
  numls = -1;
    
  # Print the column labels for reporting the algorithm's progress.
  if verbose
    @printf("iter      objective -min(g+1)  #nz #qp #ls\n")
  end

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  for i = 1:maxiter

    # Start timing the iteration.
    tic();
      
    # Compute the gradient and Hessian, optionally using the partial
    # QR decomposition to increase the speed of these computations.
    # gradient and Hessian computation -- Rank reduction method
    if (lowrank_rank > 0) & (gradient_calc == "approx")
        if lowrank_method == "qr"
          D = 1 ./(F[:Q]*(F[:R]*(P'*x)) .+ eps_tol);
          g = -P * F[:R]' * (F[:Q]'*D)/n;
          H = P * F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q])*F[:R]*P'/n + eps_tol * Diagonal(ones(k));
        elseif lowrank_method == "svd"
          D = 1 ./(F[:U]*(S*(F[:Vt]*x)) .+ eps_tol);
          g = -F[:Vt]'*(S * (F[:U]'*D))/n;
          H = (F[:V]*S*(F[:U]'*Diagonal(D.^2)*F[:U])* S*F[:Vt])/n + eps_tol * Diagonal(ones(k));
        else
          D = 1 ./(L*x .+ eps_tol);
          g = -L'*D/n;
          H = L'*Diagonal(D.^2)*L/n + eps_tol * Diagonal(ones(k));
        end
            
    elseif (lowrank_rank > 0) & (gradient_calc == "exact")
        error("We do not provide this case since it is inefficient");
            
    elseif lowrank_rank == 0
    #    D = 1 ./(L*x .+ eps_tol);
    #    g = -(L'*D)/n;
    #    H = speye(k)/n;
    #else
        error("We do not provide this case when lowrank_rank = 0");
    end
        
            
    if nullprior > 0
      g[1] -= nullprior/x[1]/n;
      H[1,1] += nullprior/x[1]^2/n;
    end

    # Report on the algorithm's progress.
    if lowrank_method == "qr"
      obj[i] = -sum(log.(F[:Q]*(F[:R]*(P'*x)) .+ eps_tol));
    elseif lowrank_method == "svd"
      obj[i] = -sum(log.(F[:U]*(S*(F[:Vt]*x)) .+ eps_tol));
    else
      obj[i] = mixobjective(L,x,eps_tol);
    end
    gmin[i] = minimum(g) .+ 1;
    nnz[i]  = length(find(x .> sparse_tol));
    nqp[i]  = j;
    if verbose
      if i == 1
          @printf("%4d %0.8e %+0.2e %4d\n",
              i,obj[i],-gmin[i],nnz[i]);
      else
          @printf("%4d %0.8e %+0.2e %4d %3d %3d\n",
              i,obj[i],-gmin[i],nnz[i],nqp[i-1],numls);
      end
    end
      
    # Check convergence of outer loop
    if gmin[i] >= -convergence_tol
      break
    end
    
    # Solve the QP subproblem using either the active-set or
    # interior-point (MOSEK) method.
    out, qptiming[i], bytes, gctime,
    memallocs = @timed if true
      y,nqp[i] = qpactiveset(x,g,H,convergence_tol = convergence_tol,sparse_tol = sparse_tol,
                      maxiter = maxqpiter);
    end
    
    # Perform backtracking line search
    if linesearch == true
        for t = 1:10
          if lowrank_method == "qr"
            D_new = 1 ./(F[:Q]*(F[:R]*(P'*y)) .+ eps_tol);
          elseif lowrank_method == "svd"
            D_new = 1 ./(F[:U]*(S*(F[:Vt]*y)) .+ eps_tol);
          else
            D_new = 1 ./(L*y .+ eps_tol);
          end
          if all(D_new .> 0)
            if sum(log.(D)) - sum(log.(D_new)) > sum((x-y) .* g) / 2
              break
            end
          end
          y = (y-x)/2 + x;
        end
        numls = t;
    else
        numls = -1;
    end

    # Update the solution to the original optimization problem.
    x = copy(y);

    # Get the elapsed time for the ith iteration.
    timing[i] = toq();
  end

  # Return: (1) the solution (after zeroing out any values below the
  # tolerance); (2) the value of the objective at each iteration; (3)
  # the minimum gradient value of the modified objective at each
  # iteration; (4) the number of nonzero entries in the vector at each
  # iteration; and (5) the number of inner iterations taken to solve
  # the QP subproblem at each outer iteration.
  x[x .< sparse_tol] .= 0;
  x             = x/sum(x);
  totaltime     = lowranktime + sum(timing[1:i]);
  if verbose
    @printf("Optimization took %d iterations and %0.4f seconds.\n",i,totaltime)
  end
  return Dict([("x",Vector(x)), ("totaltime",totaltime),
               ("lowranktime",lowranktime), ("obj",obj[1:i]),
               ("gmin",gmin[1:i]), ("nnz",nnz[1:i]),
               ("nqp",nqp[1:i]), ("timing",timing[1:i]),
               ("qptiming",qptiming[1:i])])
end

# Solve the QP subproblem for the mix-SQP algorithm using an active
# set method.
function qpactiveset(x, g, H; convergence_tol = 1e-8, sparse_tol = 1e-3, maxiter = 100)

  # Get the number of degrees of freedom in the optimization problem.
  k = length(x);
    
  # Initialize the solution to the QP subproblem.
  y      = sparse(zeros(k));
  ind    = find(x .> sparse_tol);
  y[ind] .= 1/length(ind);
  i = 0;

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  for i = 1:maxiter
        
    # Define the smaller QP subproblem.
    s   = length(ind);
    H_s = H[ind,ind];
    d   = H*y + 2*g .+ 1;
    d_s = d[ind];

    # Solve the smaller problem.
    p      = sparse(zeros(k));
    p_s    = -H_s\d_s;
    p[ind] = p_s;

    # Check convergence using KKT
    if norm(p_s) < convergence_tol
            
      # Compute the Lagrange multiplier.
      z = d;
      if all(z .>= -convergence_tol)
        break
      elseif length(ind) < k
        notind  = setdiff(1:k,ind);
        ind_min = notind[findmin(z[notind])[2]];
        ind     = sort([ind; ind_min]);
      end
    else
          
      # find a feasible step length.
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
          deleteat!(ind,find(ind .- ind_block .== 0));
        end
      end
          
      # Move to the new "inner loop" iterate (y) along the search
      # direction.
      y = y + alpha * p;
    end
  end

  # Return the solution to the quadratic program.
  return y, Int(i)
end