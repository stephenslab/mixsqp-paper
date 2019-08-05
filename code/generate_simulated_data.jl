# Generate the data sets used in the comparisons of mix-SQP vs. the
# first-order methods (EM and projected gradient descent).

# SCRIPT PARAMETERS
# -----------------
n        = 20000;
m        = 20;
filename = "simdata-n=20000-m=20.csv";

# SET UP ENVIRONMENT
# ------------------
using Distributions
include("datasim.jl");
include("likelihood.jl");

# Initialize the sequence of pseudorandom numbers.
srand(2019);

# SIMULATE DATA
# -------------
# Generate the n x m conditional likelihood matrix.
z  = normtmixdatasim(n);
sd = autoselectmixsd(z,nv = m);
L  = normlikmatrix(z,sd = sd);

# Save the matrix to a CSV file.
writecsv(filename,L);
