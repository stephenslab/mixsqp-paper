# make small data
# set n and n_half
n = 5 * 10^3;
n_half = round(Int,n/2);
srand(4)
x = [randn(n_half);3*randn(n-n_half)];
s = ones(n); 
temp = genL(x,s, mult = 1.35);
L1 = temp[:,1:20];

# make large data
# set n and n_half
n = 10^5;
n_half = round(Int,n/2);
srand(1)
x = [randn(n_half);3*randn(n-n_half)];
s = ones(n); 
temp = genL(x,s, mult = 1.04);
L2 = temp[:,1:100];
