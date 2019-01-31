function mixGD(L;  w = ones(size(L,2))/size(L,2), maxiter = 10000,
               alpha = 0.1, tol = 1e-6, eps = 1e-8, lowrank = "none",
               r = 1, formulation = "simplex")
    
    # Get the number of rows (n) and columns (k) of the conditional
    # likelihood matrix.
    n = size(L,1);
    k = size(L,2);
    
    if lowrank == "qr"
      F = pqrfact(L, rank = r);
      R = F[:R][:,sortperm(F[:p])];
    end

    # Initialize storage for outputs obj, maxd and timing.
    obj    = zeros(maxiter);
    maxd   = zeros(maxiter);
    timing = zeros(maxiter);
    
    # loop start
    for iter = 1:maxiter
        
        tic();
        
        # calculate gradient and objective value
        D         = 1 ./ (L * w + eps);
        obj[iter] = sum(log.(D));
        g         = (L'*D)/n;
        
        # get step
        if lowrank == "none"
          p = alpha * g;
        end
        
        # F[:Q] * F[:R][:,sortperm(F[:p])]
        # compute quasi-Hessian with lowrank r
        if lowrank == "qr"
            H     = R' * ((D .* F[:Q])' * (D .* F[:Q]))  * R / n + eps * speye(k);
            p     = alpha * (H\g);
        end
            
        # take initial step
        wnew = take_step(w, p; formulation = formulation);
        
        # line search
        for i = 1:10
            fnew = 1 ./ (L * wnew + eps);
            if obj[iter] - fnew[1] > dot(w - wnew, g) / 2;
                break;
            else
                p    = p/2;
                wnew = take_step(w, p);
            end
        end
        
        # save diff
        maxd[iter] = maximum(abs.(w - wnew));
        timing[iter] = toq();
        
        # convergence check
        if maxd[iter] < tol
          break;
        end
        
        # take step
        w = copy(wnew);
    end

    # Return the mixture weights and other optimization info.
  return w, obj[1:iter], maxd[1:iter], iter, timing[1:iter]
end
    
function take_step(w, g; formulation = "simplex");
    if formulation == "simplex"
        w = proj_simplex(w + g);
    elseif formulation == "nonnegative"
        w = max.(w + g - 1, 0);
    elseif formulation == "square"
        w = w .* g.^2;
        w = w/sum(w);
    else
        error("formulation not supported");
    end
    
    return w
end

function proj_simplex(x)
    u = sort(x, rev = true);
    sv = cumsum(u);
    rho = find(u .> (sv - 1) ./ (1:length(u)))[end]
    theta = max.(0.0, (sv[rho] - 1) ./ rho)
    return max.(x - theta, 0);
end
