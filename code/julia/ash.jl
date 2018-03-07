using LowRankApprox

include("mixSQP_time.jl")

function ash(x,s; mult = 1.3, lowrank = "svd")
    tic();
    s2 = s.^2
    x2 = x.^2;
    s_min = sqrt(minimum(s2))/10;
    if all(x2 .<= s2) s_max = 10 * s_min; # to deal with the occassional odd case
    else s_max = 2 * sqrt(maximum(x2-s2)); # this computes a rough largest value you'd want to use
    end
    
    # choose grid
    if mult == 0
        grid = [0;s_max];
    else 
        n = ceil(Int, log2(s_max/s_min)/log2(mult));
        grid = [0;mult.^((-n):0) * s_max];
    end
    
    # matrix likelihood
    s_matrix = sqrt.((s2) .+ (grid.^2)') # n by m matrix of standard deviations of convolutions
    log_lik = -(x./s_matrix).^2/2 - log.(s_matrix) - log(2*pi)/2;
    log_lik = log_lik - repmat(maximum(log_lik,2),1,size(log_lik,2));
    log_lik = exp.(log_lik);
    t1 = toq();
    
    # fit the model
    temp = mixSQP_time(log_lik, ptol = 1e-8, lowrank = lowrank);
    t3 = temp[3];
    t2 = temp[2];
    p = temp[1];
    
    tic();
    # exploit sparsity
    ind = find(p .> 0);
    ps2 = grid[ind].^2;
    
    # posterior calculation
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
    return post_mean, post_mean2, log_lik, p, [t1;t2;t3;t4]

end
