# Simulate n random numbers generated as follows: half of the random
# numbers are drawn from the standard univariate normal, and half of
# the random numbers are drawn from a univariate normal with zero mean
# and a standard deviation of 3. The return value is a vector of
# floats of length n.
function normdatasim(n::Int)

  # Check input "n".
  if n <= 0
    throw(ArgumentError("Argument \"n\" should be positive"))
  end

  # Generate the random numbers.
  n1 = round(Int,0.5*n);
  n2 = round(Int,0.2*n);
  n3 = n - n1 - n2;
  return vcat(randn(n1),4*randn(n2),6*randn(n3))
end

# Simulate n random numbers generated as follows: 50% of the random
# numbers are drawn from the standard univariaten normal; 20% are
# drawn from a t-distribution with 4 degrees of freedom; and the
# remaining 30% are drawn from a t-distribution with 6 degrees of
# freedom. The return value is a vector of floats of length n.
function normtmixdatasim(n::Int)
  
  # Check input "n".
  if n <= 0
    throw(ArgumentError("Argument \"n\" should be positive"))
  end

  # Generate the random numbers.
  n1 = round(Int,0.5*n);
  n2 = round(Int,0.2*n);
  n3 = n - n1 - n2;
  return vcat(randn(n1),rand(TDist(4),n2),rand(TDist(6),n3))
end
