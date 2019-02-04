% TO DO: Explain here what this script does.
rng(1);

% Simulate a data set.
n = 1000;
p = 100;
A = randn(n,p);
x = rand(p,1) .* (rand(p,1) > 0.5);
b = A*x + randn(n,1);

% Fit a linear regression subject to the coefficients being non-negative.
lb = zeros(p,1);
ub = Inf(p,1);
x0 = zeros(p,1);
f  = @(x) SquaredError(x,A,b);
g  = @(x) boundProject(x,lb,ub);
out = minConf_SPG(f,x0,g);

