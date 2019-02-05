% Small illustration of the projected gradient algorithm applied to the
% problem of fitting a linear regression subject to the coefficients being
% probabilities (that is, they are non-negative, and sum to one). This
% example is drawn from the code on Mark Schmidt's website.
rng(1);

% SIMULATE DATA
% -------------
n = 1e4;
p = 10;
X = randn(n,p);
w = rand(p,1) .* (rand(p,1) > 0.75);
w = w/sum(w);
y = X*w + randn(n,1);

% FIT MODEL
% ---------
w0 = zeros(p,1);
f  = @(w) SquaredError(w,X,y);
g  = @(w) projectSimplex(w);
wspg = minConf_SPG(f,w0,g,struct('useSpectral',true,'optTol',1e-8,...
                                 'progTol',1e-15,'suffDec',1e-8));
[w wspg]