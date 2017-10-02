include("sqp.jl")

L = readdlm("normlik.5000x20.txt",' ');
k = size(L,2)

x   = ones(k)/k;
out = sqp(L,x,verbose = true)
