# Try to select a reasonable set of sigma values that should be used
# for the adaptive shrinkage model based on the values of x (the noisy
# observations) and s (the standard error in the observations). The
# return value is a vector of sigma values.
function autoselectmixsd(x::Array{Float64,1},
                         s::Array{Float64,1} = ones(size(x));
                         gridmult::Float64 = 1.4,
                         nv::Int = 0)

  # Get the number of samples.
  n = length(x);
    
  # Check input "s"---it should be the same length.
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

  # Check input "nv".
  if !(nv == 0 || nv > 1)
    throw(ArgumentError("Input \"nv\" should be 0, or greater than 1"))
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
  if nv > 0
    return vcat(0,logspace(log10(smin),log10(smax),nv - 1));
  elseif gridmult == 0
    return vcat(0,smax)
  else 
    m = ceil(Int,log2(smax/smin)/log2(gridmult));
    return vcat(0,logspace(log10(smin),log10(smax),m + 1));
  end
end

# Compute the n x k conditional likelihood matrix, where n is the
# number of samples and k is the number of mixture components, for the
# case when the likelihood is univariate normal and prior is a mixture
# of univariate normals.
#
# Entry (i,j) of the conditional likelihood matrix is equal to
# N(0,s[i]^2 + sd[j]^2), the normal density with zero mean and
# variance s[i]^2 + sd[j]^2.
#
# If normalize = true, each row of the likelihood matrix is divided by
# the largest value in the row. After normalization, the largest value
# in each row is 1.
function normlikmatrix(x::Array{Float64,1},
                       s::Array{Float64,1} = ones(size(x));
                       sd::Array{Float64,1} = autoselectmixsd(x,s),
                       normalizerows::Bool = true)

  # Get the number of samples (n) and the number of prior mixture
  # components (k).
  n = length(x);
  k = length(sd);
    
  # Check input "s"---it should be the same length at x.
  if length(s) != n
    throw(ArgumentError("Arguments \"x\" and \"s\" should have the same" *
                        "length"))
  elseif any(s .<= 0)
    throw(ArgumentError("All elements of \"s\" should be positive"))
  end

  # Check input "sd".
  if any(sd .< 0)
    throw(ArgumentError("All elements of \"sd\" should be non-negative"))
  end
    
  # Compute the n x k matrix of standard deviations.
  S = sqrt.((s.^2) .+ (sd.^2)');

  # Compute the log-densities, and normalize the rows, if requested.
  L = -(x./S).^2/2 - log.(S) .- log(2*pi)/2;
  if normalizerows

    # This is the same as
    #
    #   L = L - repmat(maximum(L,2),1,k);
    #
    # but uses memory more efficiently to complete the operation.
    L = broadcast(-,L,maximum(L,dims = 2));
  end
  return exp.(L)
end

# Implements the "logspace" function that was available in an
# earlier version of Julia.
function logspace(a, b, n)
  return exp10.(range(a,stop = b,length = n))
end

