using LowRankApprox

function sqp_conv_analysis(L,eps,tol,sptol,maxiter = 100, option = true)
  n = size(L,1); k = size(L,2);
  F = pqrfact(L, rtol=tol);
  P = sparse(F[:P]);
  iter = maxiter;
  # initialize
  x = zeros(k); x[1] = 1/2; x[k] = 1/2; D = 0; i = 0;
  f = zeros(iter);
  l = zeros(iter);

  # QP subproblem start
  for i = 1:iter
    # gradient and Hessian computation -- Rank reduction method
    D = 1./(F[:Q]*(F[:R]*(P'*x)) + eps);
    g = -P * F[:R]' * (F[:Q]'*D)/n;
    H = P * F[:R]' * (F[:Q]'*Diagonal(D.^2)*F[:Q]) * F[:R] * P'/n + eps * eye(k);
    f[i] = sum(log(max(D,eps)))/n;
    # initialize
    ind = find(x .> sptol);
    y = sparse(zeros(k)); y[ind] = 1/length(ind);
    l[i] = sum((g[ind] + 1) .* x[ind]);

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
    if minimum(g+1) >= 0 && option == false
      break;
    end
  end
  f = f[2:i];
  l = l[2:i];
  x[x .< sptol] = 0;
  return full(x), f, l
end
