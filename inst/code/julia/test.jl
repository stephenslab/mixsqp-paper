using Plots
include("mixsqp.jl")

# Load the conditional likelihood matrix.
L = readdlm("../../datafiles/normlik.5000x20.txt",' ');
k = size(L,2)

# Choose an initial solution.
x = ones(k)/k;

# Run the SQP algorithm.
out = mixsqp(L,x);

# The the SQP algorithm using the partial QR decmoposition.
out2 = mixsqp(L,x,pqrtol = 1e-8);

# Show the progress toward the minimum.
n = length(out[2])
plot(1:n,log.(out[2] - minimum(out[2])),linewidth = 2)
