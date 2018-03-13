# TO DO: Edit comments for this function.
#
# try to select a default range for the sigmaa values that should be
# used, based on the values of betahat and sebetahat mode is the
# location about which inference is going to be centered mult is the
# multiplier by which the sds differ across the grid grange is the
# user-specified range of mixsd.
#
function autoselectmixsd(x::Array{Float64,1},
                         s::Array{Float64,1} = ones(size(x)),
                         gridmult::Float64 = 1.4)

  # Check input "s".
  # TO DO.
    
  # Check input "gridmult".
  # TO DO.
  x2 = x.^2;
  s_min = sqrt(minimum(s.^2))/10;
    if all(x2 .<= s2)

      # This deals with the occasional boundary case.
      smax = 10 * smin; 
    else

      # This is, roughly, the largest value you'd want to use.
      s_max = 2 * sqrt(maximum(x2-s2)); 
    end
    
  # Choose the grid of sigmas.
  if gridmult == 0
    grid = vcat(0,smax)
  else 
        n = ceil(Int, log2(s_max/s_min)/log2(mult));
        grid = [0;mult.^((-n):0) * s_max];
    end

end


function normlikmatrix (x::Array{Float64,1},
                        s::Array{Float64,1} = ones(size(x)),
                        σ::)
    s2 = s.^2
    x2 = x.^2;
    
    # matrix likelihood
    s_matrix = sqrt.((s2) .+ (grid.^2)') # n by m matrix of standard deviations of convolutions
    L = -(x./s_matrix).^2/2 - log.(s_matrix) - log(2*pi)/2;
    L = L - repmat(maximum(L,2),1,size(L,2));
    L = exp.(L);
    
    return L
end
