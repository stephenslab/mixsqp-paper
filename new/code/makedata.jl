using DataFrames

# load ash
include("ash.jl")

# make large data
srand(1)
x = [randn(50000);3*randn(50000)];
s = ones(10^5); 
temp = ash(x,s, mult = 1.04);
L2 = temp[3][:,1:100];

# save
writetable("data/sample100000x100.txt", DataFrame(L2), header = false, separator = ' ');