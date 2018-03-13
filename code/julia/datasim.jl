# Simulate n random numbers generated as follows: half of the random
# numbers are drawn from the standard univariate normal, and half of
# the random numbers are drawn from a univariate normal with zero mean
# and a standard deviation of 3. The return value is a vector of
# floats of length n.
function normdatasim(n::Int)

  # Check input "n".
  if n <= 0
    throw(ArgumentError("Argument \"n\" must be a positive integer"))
  end

  # Generate the random numbers.
  n1 = round(Int,n/2);
  n2 = n - n1;
  return vcat(randn(n1),3*randn(n2))
end
