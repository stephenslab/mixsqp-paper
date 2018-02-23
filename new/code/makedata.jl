using DataFrames

# load ash
include("ash.jl")

# set n and n_half
n = 10^5;
n_half = round(Int,n/2);

# make large data
srand(1)
x = [randn(n_half);3*randn(n-n_half)];
s = ones(n); 
temp = ash(x,s, mult = 1.04);
L2 = temp[3][:,1:100];

# save
writetable("data/sample100000x100.txt", DataFrame(L2), header = false, separator = ' ');