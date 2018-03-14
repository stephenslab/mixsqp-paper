function ash(x,s; method = "mixSQP", mult = 1.3, lowrank = "svd")

    # Compute the n x m likelihood matrix.
    tic();
    sd      = autoselectmixsd(x,s,gridmult = mult);
    log_lik = normlikmatrix(x,s,sd = sd);
    t1      = toq();
    
    # Fit the model.
    if method == "mixSQP"
        temp = mixSQP_time(log_lik, pqrtol = 1e-8, lowrank = lowrank);
        t3 = temp[3];
        t2 = temp[2];
        p = temp[1];
    elseif method == "REBayes"
        temp = REBayes(log_lik);
        t2 = temp[2];
        p = temp[1];
    else
        error("Error")
    end
    
    tic();
    
    # Exploit sparsity.
    ind = find(p .> 0);
    ps2 = grid[ind].^2;
    
    # Posterior calculations.
    temp = s2 .+ ps2';
    comp_post_mean = (x * ps2') ./ temp;
    comp_post_sd2 = (s2 * ps2') ./ temp;
    comp_post_mean2 = comp_post_sd2 + comp_post_mean.^2;
    comp_post_prob = log_lik[:,ind] .* p[ind]';
    comp_post_prob = comp_post_prob ./ sum(comp_post_prob,2);
    post_mean = sum(comp_post_prob .* comp_post_mean,2);
    post_mean2 = sum(comp_post_prob .* comp_post_mean2,2);
    t4 = toq();
    
    # return posterior first/second moments
    if method == "mixSQP"
        return post_mean, post_mean2, log_lik, p, [t1;t2;t3;t4]
    else
        return post_mean, post_mean2, log_lik, p, [t1;t2;t4]
    end

end
