# TO DO: Edit comments for this function.
#
# try to select a default range for the sigmaa values that should be
# used, based on the values of betahat and sebetahat mode is the
# location about which inference is going to be centered mult is the
# multiplier by which the sds differ across the grid grange is the
# user-specified range of mixsd.
#
function autoselectmixsd(x::Array{Float64,1},
                         s::Array{Float64,1} = ones(size(x));
                         gridmult::Float64 = 1.4)

  # Get the number of samples.
  n = length(x);
    
  # Check input "s"---it should be the same length
  if length(s) != n
    throw(ArgumentError("Arguments \"x\" and \"s\" should have the same" *
                        "length"))
  elseif any(s .<= 0)
    throw(ArgumentError("All elements of \"s\" should be positive"))
  end
      
  # Check input "gridmult".
  if gridmult < 0
    throw(ArgumentError("Input \"gridmult\" should be non-negative"))
  end

  # Determine the minimum and maximum sigma settings.
  smin = sqrt(minimum(s.^2))/10;
  if all(x.^2 .<= s.^2)

    # This deals with the occasional boundary case.
    smax = 10 * smin; 
  else

    # This is, roughly, the largest value you'd want to use.
    smax = 2 * sqrt(maximum(x.^2 - s.^2)); 
  end
    
  # Choose the grid of sigmas.
  if gridmult == 0
    return vcat(0,smax)
  else 
    m = ceil(Int,log2(smax/smin)/log2(gridmult));
    return vcat(0,gridmult.^colon(-m,0) * smax)
  end
end

## function normlikmatrix(x::Array{Float64,1},
##                        s::Array{Float64,1} = ones(size(x)),
##                        sigma)
##     s2 = s.^2
##     x2 = x.^2;
    
##     # matrix likelihood
##     s_matrix = sqrt.((s2) .+ (grid.^2)') # n by m matrix of standard deviations of convolutions
##     L = -(x./s_matrix).^2/2 - log.(s_matrix) - log(2*pi)/2;
##     L = L - repmat(maximum(L,2),1,size(L,2));
##     L = exp.(L);
    
##     return L
## end
