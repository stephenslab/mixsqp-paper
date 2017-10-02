using Plots
include("sqp.jl")

# Load the conditional likelihood matrix.
L = readdlm("normlik.5000x20.txt",' ');
k = size(L,2)

# Choose an initial solution.
x = ones(k)/k;

# Run the SQP algorithm.
out = sqp(L,x,verbose = true)

# Show the progress toward the minimum.
n = length(out[2])
plot(1:n,log.(out[2] - minimum(out[2])),linewidth = 2)
