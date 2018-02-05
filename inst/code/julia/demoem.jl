# TO DO: Explain here what this script does, and how to use it.
using HDF5
using Plots
include("mixem.jl")

# Load the 5000 x 20 conditional likelhood matrix computed from a
# simulated data set.
filename = "../../datafiles/normmix.data.h5";
L   = h5read(filename,"L");
w   = h5read(filename,"w");
obj = h5read(filename,"obj")[1];

# Fit the model using the EM algorithm.
k   = size(L,2);
w0  = ones(k)/k;
wem, objem, maxd = mixem(L,w0,tol = 1e-5);

# This plot shows that convergence to the minimum is very slow.
plot(1:length(objem),log10(objem - obj))
