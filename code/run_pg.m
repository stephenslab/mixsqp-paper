% This script completes the comparison of mix-SQP against the first-order
% methods; compare_first_order.jl runs EM and the mix-SQP variants, and
% this script runs the projected gradient method.
% 
% To run on the RCC cluster, I set up my environment with the
% following commands:
%
%   sinteractive --partition=broadwl --mem=8G
%   module load matlab/2018b
%

% Load the data.
fprintf('Reading data.\n')
L = csvread('simdata-n=2000-m=200.csv');
[n m] = size(L);

% Fit the model using the projected gradient algorithm.
fprintf('Fitting model using projected gradient method.\n');
x0 = ones(m,1)/m;
fx = @(x) mixobj(L,x,1e-8);
g  = @(x) projectSimplex(x);
[x f nfevals nproj timings] = ...
    minConf_SPG(fx,x0,g,struct('useSpectral',false,'optTol',1e-8,...
                               'progTol',1e-15,'suffDec',1e-8,...
                               'memory',1,'maxIter',1e4,'verbose',0,...
                               'interp',2,'testOpt',1,'curvilinear',0));

% Write the objective value and runtime at each iteration to a CSV file.
fprintf('Writing results to file.\n');
csvwrite('pg-n=2000-m=20.csv',[f/n timings]);