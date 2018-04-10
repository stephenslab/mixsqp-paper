# TO DO: Explain here what this script is for.
using Distributions
using LowRankApprox
include("datasim.jl")
include("likelihood.jl")

# Initialize the sequence of pseudorandom numbers.
srand(1);

# Generate a large data set.
x = normtmixdatasim(round(Int,1e5));

# Compute the n x k likelihood matrix for a mixture of zero-centered
# normals, with k = 100.
@printf "Determine variance parameters for mixture components.\n"
@time sd = autoselectmixsd(x,nv = 100);
@printf "Compute likelihood matrix, L.\n"
@time log_lik = normlikmatrix(x,sd = sd);

# Compute partial QR factorization of the likelihood matrix.
@printf "Compute partial QR factorization of L.\n"
@time F = pqrfact(L,rtol = 1e-12);

# Here we get a more detailed breakdown of the likelihood computations.
#
# Compute the n x k matrix of standard deviations.
@printf "Breakdown of likelihood computations:\n"
s = ones(size(x));
@printf "S = sqrt.((s.^2) .+ (sd.^2)') \n"
@time S = sqrt.((s.^2) .+ (sd.^2)');
@printf "L = -(x./S).^2/2 \n"
@time L = -(x./S).^2/2;
@printf "L = L - log.(S) \n"
@time L = L - log.(S)
@printf "L = L - log(2*pi)/2 \n"
@time L = L - log(2*pi)/2;
@printf "L = broadcast(-,L,maximum(L,2)) \n"
@time L = broadcast(-,L,maximum(L,2));
