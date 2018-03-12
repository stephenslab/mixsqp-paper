function genL(x,s; mult = 1.4)
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
    
    return log_lik
end
