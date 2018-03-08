using Mosek, JuMP

function QPsubprob(L; initialize = "equal", method = "activeset", warmstart = true, lowrank = "svd")
  # fix setting
  eps = 1e-8; tol = 1e-8; pqrtol = 1e-10; sptol=1e-3; maxiter = 100;
  n = size(L,1); k = size(L,2);
  if lowrank == "qr"
      F = pqrfact(L, rtol=pqrtol);
      P = sparse(F[:P]);
  elseif lowrank == "svd"
      F = psvdfact(L, rtol=pqrtol);
      S = Diagonal(F[:S]);
  else
  end
    
  # return param
  timing = zeros(maxiter); linesearch = zeros(maxiter); i = 0;
    
  # initialize
  if initialize == "equal"
      x = ones(k)/k;
  elseif initialize == "null"
      x = sparse(zeros(k)); x[1] = 1;
  elseif initialize == "sparse"
      x = sparse(zeros(k)); x[1] = 1/2; x[round(Int,k/2)] = 1/2;
  else
      x = rand(k); x = x/sum(x);
  end

  # QP subproblem start
  for i = 1:maxiter
    # gradient and Hessian computation -- Rank reduction method
    if lowrank == "qr"
        D = 1./(F[:Q]*(F[:R]*(P'*x)) + eps);
        g = -P * F[:R]' * (F[:Q]'*D)/n;
        H = P * F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q]) * F[:R] * P'/n + eps * eye(k);
    elseif lowrank == "svd"
        D = 1./(F[:U]*(S*(F[:Vt]*x)) + eps);
        g = -F[:Vt]'*(S * (F[:U]'*D))/n;
        H = (F[:V]*S*(F[:U]'*Diagonal(D.^2)*F[:U])* S*F[:Vt])/n + eps * eye(k);
    else
        D = 1./(L*x + eps);
        g = -L'*D/n;
        H = L'*Diagonal(D.^2)*L/n + eps * eye(k);
    end
    
    tic();
    if method == "activeset"
        # initialize working set warmstart or not
        if warmstart == true
            ind = find(x .> sptol); # this is inactive set
            y = sparse(zeros(k)); y[ind] = 1/length(ind);
        else
            ind = [1];
            y = sparse(zeros(k)); y[1] = 1;
        end
            
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
    else
        # interior point method start
        mod = Model(solver=MosekSolver(QUIET = true));
        @variable(mod, y[1:k] >= 0);
        @NLobjective(mod, Min, sum(y[i]'*H[i,j]*y[j] for i = 1:k, j = 1:k)/2 + sum(y[i]*(2*g[i]+1) for i = 1:k));
        solve(mod);
        y = getvalue(y);
    end
    timing[i] = toq();
    
     
    for linesearch[i] = 1:10
        if sum(log.(D)) - sum(log.(1./(F[:U]*(S*(F[:Vt]*y)) + eps))) > sum((x-y) .* g) / 2
            break;
        end
        y = (y-x)/2 + x;
    end
    x = y;

    # convergence check
    if minimum(g+1) >= 0
      break;
    end
  end
  x[x .< sptol] = 0; x = x/sum(x);
    
  return Dict([("numiter", i), ("eachQPtime",timing[1:i]), ("linesearch", linesearch[1:i]),
            ("totalQPtime", sum(timing[1:i]))])
end