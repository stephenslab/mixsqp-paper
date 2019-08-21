# Outputs the value of the objective function at x.
function mixobjective(L, x, eps = 0)
  return -sum(log.(L * x .+ eps))
end
          
# Fit a mixture model using the sequential quadratic programming
# algorithm ("mix-SQP"). Input argument L is the n x m conditional
# likelihood matrix, where n is the number of samples and m is the
# number of mixture components. Optional input argument x is the
# initial estimate of the mixture weights (when one or more entries of
# x are negative, the default initial estimate is used in which all
# the entries are set to 1/k, where k = size(L,2).
function mixSQP(L; x = -1, convtol = 1e-8, pqrtol = 1e-8, eps = 1e-8,
                sptol = 1e-6, maxiter = 100, maxqpiter = 100,
                lowrank = "svd", qpsubprob = "activeset",
                verbose = true)
    
  # Get the number of rows (n) and columns (k) of L.
  n = size(L,1);
  k = size(L,2);

  # If any entries of input x are negative, set the initial estimate
  # to the default setting. When the MOSEK solver is used, the initial
  # estimate is set to a sparse vector in which all but two of the
  # entries are zero.
  if any(x .< 0)
    x = ones(k)/k;
  end
    
  # If requested (i.e., if pqrtol > 0), compute the partial QR
  # decomposition with relative precision "tol", then retrieve the
  # permutation matrix, P. For details on the partial QR, see
  # https://github.com/klho/LowRankApprox.jl.
  if lowrank != "none" && qpsubprob == "mosek"
    error("lowrank must be \"none\" when qpsubprob = \"mosek\"");
  end
  if lowrank == "qr"
    F = pqrfact(L, rtol = pqrtol);
    P = sparse(F[:P]);
  elseif lowrank == "svd"
    F = psvdfact(L, rtol = pqrtol);
    S = Diagonal(F[:S]);
  end
    
  # Summarize the analysis here.
  if verbose
    @printf("Running SQP algorithm with the following settings:\n")
    @printf("- %d x %d data matrix\n",n,k)
    @printf("- convergence tolerance  = %0.2e\n",convtol)
    @printf("- zero threshold         = %0.2e\n",sptol)
    if lowrank == "qr"
      err = maximum(abs.(F[:Q]*F[:R]*P' - L));
      @printf("- partial QR tolerance   = %0.2e\n",pqrtol)
      @printf("- partial QR max. error  = %0.2e\n",err)
    elseif lowrank == "svd"
      err = maximum(abs.(F[:U]*S*F[:Vt] - L));
      @printf("- partial SVD tolerance  = %0.2e\n",pqrtol)
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
    
  # Print the column labels for reporting the algorithm's progress.
  if verbose
    @printf("iter      objective -min(g+1)  #nz #qp #ls\n")
  end

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  numiter = 0;
  numls   = 0;
  for i = 1:maxiter
    numiter = i;
      
    out, timing[i] = @timed begin
      
      # Compute the gradient and Hessian, optionally using the partial
      # QR decomposition to increase the speed of these computations.
      # gradient and Hessian computation -- Rank reduction method
      if lowrank == "qr"
        D = 1 ./ (F[:Q]*(F[:R]*(P'*x)) .+ eps);
        g = -P * (F[:R]' * (F[:Q]'*D)/n);
        H = P * (F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q])*F[:R])*P'/n + 
            eps*Diagonal(ones(k));
      elseif lowrank == "svd"
        D = 1 ./ (F[:U]*(S*(F[:Vt]*x)) .+ eps);
        g = -F[:Vt]'*(S * (F[:U]'*D))/n;
        H = (F[:V]*S*(F[:U]'*Diagonal(D.^2)*F[:U])* S*F[:Vt])/n + 
            eps*Diagonal(ones(k));
      else
        D = 1 ./ (L*x .+ eps);
        g = -L'*D/n;
        H = L'*Diagonal(D.^2)*L/n .+ eps * Diagonal(ones(k));
      end
         
      # Report on the algorithm's progress.
      obj[i]  = mixobjective(L,x,eps);
      gmin[i] = minimum(g) .+ 1;
      nnz[i]  = length(findall(x .> sptol));
      if verbose
        if i == 1
          @printf("%4d %0.8e %+0.2e %4d\n",i,obj[i],-gmin[i],nnz[i]);
        else
          @printf("%4d %0.8e %+0.2e %4d %3d %3d\n",i,obj[i],
                  -gmin[i],nnz[i],nqp[i - 1],numls);
        end
      end
      
      # Check convergence of outer loop.
      if gmin[i] >= -convtol
        break
      end
    
      # Solve the QP subproblem using either the active-set or
      # interior-point (MOSEK) method.
      out, qptiming[i] = @timed if qpsubprob == "activeset"
        y, nqp[i] = qpactiveset(x,g,H,convtol = convtol,sptol = sptol,
                                maxiter = maxqpiter);
      elseif qpsubprob == "mosek"
        y = qpmosek(x,g,H);
      end
    
      # Perform backtracking line search
      for t = 1:10
        numls = t;
        if lowrank == "qr"
          D_new = 1 ./ (F[:Q]*(F[:R]*(P'*y)) .+ eps);
        elseif lowrank == "svd"
          D_new = 1 ./ (F[:U]*(S*(F[:Vt]*y)) .+ eps);
        else
          D_new = 1 ./ (L*y .+ eps);
        end
        if all(D_new .> 0)
          if sum(log.(D)) - sum(log.(D_new)) > sum((x - y) .* g) / 2
            break
          end
        end
        y = (y-x)/2 + x;
      end

      numiter = i;
    end
      
    # Update the solution to the original optimization problem.
    x = copy(y);
  end

  # Return: (1) the solution (after zeroing out any values below the
  # tolerance); (2) the value of the objective at each iteration; (3)
  # the minimum gradient value of the modified objective at each
  # iteration; (4) the number of nonzero entries in the vector at each
  # iteration; (5) the number of inner iterations taken to solve the
  # QP subproblem at each outer iteration; (6) the total runtime for
  # each iteration; and (7) the runtime for solving the quadratic
  # subproblem at each iteration.
  x[x .< sptol] .= 0;
  x = x/sum(x);
  return Dict([("x",Vector(x)),
               ("obj",obj[1:numiter]),
               ("gmin",gmin[1:numiter]),
               ("nnz",nnz[1:numiter]),
               ("nqp",nqp[1:numiter]),
               ("timing",timing[1:numiter]),
               ("qptiming",qptiming[1:numiter])])
end

# Solve the QP subproblem for the mix-SQP algorithm using an active
# set method.
function qpactiveset(x, g, H; convtol = 1e-8, sptol = 1e-6, maxiter = 100)

  # Get the number of degrees of freedom in the optimization problem.
  k = length(x);
    
  # Initialize the solution to the QP subproblem.
  y       = zeros(k);
  ind     = findall(x .> sptol);
  y[ind] .= 1/length(ind);

  # Repeat until we reach the maximum number of iterations, or until
  # convergence is reached.
  numiter = 0;
  for i = 1:maxiter
    numiter = i;
      
    # Define the smaller QP subproblem.
    s   = length(ind);
    H_s = H[ind,ind];
    d   = H*y + 2*g .+ 1;
    d_s = d[ind];

    # Solve the smaller problem.
    p      = zeros(k);
    p_s    = -H_s \ d_s;
    p[ind] = p_s;

    # Check convergence using KKT
    if norm(p_s) < convtol
            
      # Compute the Lagrange multiplier.
      z = d;
      if all(z .>= -convtol)
        break
      elseif length(ind) < k
        notind  = setdiff(1:k,ind);
        ind_min = notind[findmin(z[notind])[2]];
        ind     = sort([ind; ind_min]);
      end
    else
          
      # Find a feasible step length.
      alpha     = 1;
      alpha0    = -y[ind] ./ p_s;
      ind_block = findall(p_s .< 0);
      alpha0    = alpha0[ind_block];
      if ~isempty(ind_block)
        v, t = findmin(alpha0);
        if v < 1

          # Blocking constraint.
          ind_block = ind[ind_block[t]]; 
          alpha     = v;
              
          # Update working set if there is a blocking constraint.
          deleteat!(ind,findall(ind .- ind_block .== 0));
        end
      end
          
      # Move to the new "inner loop" iterate (y) along the search
      # direction.
      y = y + alpha * p;
    end
  end

  # Return the solution to the quadratic program.
  return y, numiter
end

# Solve the QP subproblem for the mix-SQP algorithm using MOSEK.
function qpmosek(x, g, H)
  k      = length(g);
  y      = copy(x);
  bkx    = repeat([MSK_BK_LO],k,1)[:];
  blx    = zeros(k);
  bux    = Inf * ones(k);
  numvar = length(bkx);
  c      = 2*g .+ 1;
    
  maketask() do task
    appendvars(task,numvar)
    putclist(task,[1:numvar;],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    Q    = LowerTriangular(H);
    ind  = findall(Q .> 0);
    a    = (i->i[1]).(ind)
    b    = (i->i[2]).(ind)
    putqobj(task,a,b,Q[ind])
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task)
    y = getxx(task,MSK_SOL_ITR)
  end
  return y
end
