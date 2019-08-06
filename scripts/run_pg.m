% This script completes the comparison of mix-SQP against the first-order
% methods; compare_first_order.jl runs EM and the mix-SQP variants, and
% this script runs the projected gradient method.

% SCRIPT PARAMETERS
% -----------------
datfile = 'simdata-n=20000-m=20.csv';
outfile = 'pg-n=20000-m=20.csv';

% LOAD DATA
% ---------
fprintf('Reading data.\n')
L = csvread(datfile);
[n m] = size(L);

% FIT MODEL
% ---------
% Fit the model using the projected gradient algorithm.
fprintf('Fitting model to %d x %d data set using PG method.\n',n,m);
x0 = ones(m,1)/m;
fx = @(x) mixobj(L,x,1e-8);
g  = @(x) projectSimplex(x);
[x f nfevals nproj timings] = ...
    minConf_SPG(fx,x0,g,struct('useSpectral',false,'optTol',1e-8,...
                               'progTol',1e-15,'suffDec',1e-8,...
                               'memory',1,'maxIter',5e4,'verbose',1,...
                               'interp',2,'testOpt',1,'curvilinear',0));

% Write the objective value and runtime at each iteration to a CSV file.
fprintf('Writing results to file.\n');
csvwrite(outfile,[f/n timings]);