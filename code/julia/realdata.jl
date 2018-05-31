using LowRankApprox
include("datasim.jl")
include("likelihood.jl")
for i = 1:10
    x = normdatasim(round(Int,2^i*2e4));
    sd = autoselectmixsd(x,nv = 100);
    L = normlikmatrix(x,sd = sd);
    mixSQP_record_runtime(L, tol = 1e-6)
end