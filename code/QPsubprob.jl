function QPsubprob(L; method = "activeset", eps = 1e-8, sptol = 1e-6,
                   maxiter = 1000,verbose = true)
  tol = 1e-8; 
  n = size(L,1);
  k = size(L,2);
    
  # return param
  timing = zeros(maxiter);
  linesearch = zeros(maxiter);
  y_nnz = zeros(maxiter);
  q_nnz = zeros(maxiter);
  i = 0;
    
  # initialize
  x = sparse(zeros(k));
  x[1:2] = 1/2;

  # This initialization doesn't work with MOSEK:
  # 
  #   x = sparse(ones(k)/k);
  #
    
  # Print the column labels for reporting the algorithm's progress.
  if verbose
    @printf("iter      objective -min(g+1) #nnz\n")
  end
    
  # QP subproblem start
  for i = 1:maxiter
    # gradient and Hessian computation -- Rank reduction method
    D = 1./(L*x + eps);
    g = -L'*D/n;
    H = L'*Diagonal(D.^2)*L/n + eps * eye(k);
    
    tic();
    if method == "activeset"
        y = sparse(zeros(k)); y[1] = 1; ind = [1];
            
        # Active set method start
        for j = 1:100
          # define smaller problem
          s = length(ind);
          H_s = H[ind,ind];
          d = H * y + 2 * g + 1;
          d_s = d[ind];

          # solve smaller problem
          p = sparse(zeros(k));
          p_s = -H_s\d_s; p[ind] = p_s;

          # convergence check
          if norm(p_s) < tol
            # compute the Lagrange multiplier
            lambda = d - minimum(d_s);
            # convergence test
            if minimum(lambda) >= 0;
              break;
            else
              ind_min = findmin(lambda)[2];
              ind = sort([ind;ind_min]);
            end

          # do update otherwise
          else
            # retain feasibility
            alpha = 1;
            alpha_temp = -y[ind]./p_s;
            ind_block = find(p_s .< 0);
            alpha_temp = alpha_temp[ind_block];
            if ~isempty(ind_block)
              temp = findmin(alpha_temp);
              if temp[1] < 1
                ind_block = ind[ind_block[temp[2]]]; # blocking constraint
                alpha = temp[1];
                # update working set -- if there is a blocking constraint
                deleteat!(ind, find(ind - ind_block .== 0));
              end
            end
            # update
            y = y + alpha * p;
          end
        end
    elseif method == "mosek2"
        bkx   = repmat([MSK_BK_LO], k, 1)[:]
        blx   = zeros(k)
        bux   = Inf * ones(k)
        numvar = length(bkx)
        c     = 2*g+1

        maketask() do task
            appendvars(task,numvar)
            putclist(task,[1:numvar;],c)
            putvarboundslice(task,1,numvar+1,bkx,blx,bux)
            Q = LowerTriangular(H); ind = (Q .> 0); a,b = findn(ind);
            putqobj(task,a,b,Q[ind])
            putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

            # Optimize
            optimize(task)

            y = getxx(task,MSK_SOL_ITR)
        end
   else
        mod = Model(solver=MosekSolver(QUIET = true));
        a,b = findn(H);
        @variable(mod, u[1:k] >= 0);
        @objective(mod, Min, QuadExpr(u[a],u[b],H[:]/2,AffExpr(u, 2g+1, 0)) )
        solve(mod)
        y = getvalue(u);
    end

    timing[i] = toq();
    for linesearch[i] = 1:10
        if sum(log.(D)) - sum(log.(1./(L*y + eps))) > sum((x-y) .* g) / 2
            break;
        end
        y = (y-x)/2 + x;
    end
      
    y_nnz[i] = sum(y .> sptol);
    q_nnz[i] = sum(abs.(x-y) .> sptol);
        
    x = y + 0.0;

    # convergence check
    if verbose
      obj = mixobjective(L,x,eps);
      @printf("%4d %0.8e %+0.2e %4d\n",i,obj,-minimum(g+1),sum(x .> sptol));
    end
    if minimum(g+1) >= 0
      break;
    end
  end

  x[x .< sptol] = 0; x = x/sum(x);
    
  return Dict([("x", full(x)), ("numiter", i), ("eachQPtime",timing[1:i]), ("linesearch", linesearch[1:i]),
            ("totalQPtime", sum(timing[1:i])), ("y_nnz", y_nnz[1:i]), ("q_nnz", q_nnz[1:i])])
end
