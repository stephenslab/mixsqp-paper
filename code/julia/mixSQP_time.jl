function mixSQP_time(L; eps=1e-8, tol=1e-8, ptol = 1e-10, sptol=1e-3, lowrank = "svd")
  n = size(L,1); k = size(L,2);
  tic();
  if lowrank == "qr"
      F = pqrfact(L, rtol=ptol);
      P = sparse(F[:P]);
  elseif lowrank == "svd"
      F = psvdfact(L, rtol=ptol);
      S = Diagonal(F[:S]);
  else
  end
  t = toq();
    
  iter = 100;
  # initialize
  x = zeros(k); x[1] = 1/2; x[k] = 1/2;

  tic();
  # QP subproblem start
  for i = 1:iter
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
        
    # initialize
    ind = find(x .> sptol);
    y = sparse(zeros(k)); y[ind] = 1/length(ind);

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
    x = y;

    # convergence check
    if minimum(g+1) >= 0
      break;
    end
  end
  x[x .< sptol] = 0;
  t2 = toq();
    
  return full(x/sum(x)), t, t2
end

function eval_f(L,x; eps = 1e-8)
    return -sum(log.(L*x + eps))/size(L,1) + sum(x) - 1
end

function rel_error(L,x1,x2)
    return min(abs(1-eval_f(L,x1)/eval_f(L,x2)), abs(1-eval_f(L,x2)/eval_f(L,x1)))
end
