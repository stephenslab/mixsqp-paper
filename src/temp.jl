include("sqp.jl")

L = readdlm("normlik.5000x20.txt",' ');
k = size(L,2)

x    = zeros(k);
x[1] = 0.5;
x[k] = 0.5;

out = sqp(L,x,eps = 1e-8,tol = 1e-8,sptol = 1e-8);
