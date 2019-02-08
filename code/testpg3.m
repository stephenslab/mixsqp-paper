% This is another short script to test the projected gradient method. Here
% we apply it to a less trivial maximum-likelihood estimation problem
% with a 100,000 x 10 matrix.
%
% This is the R code I used to generate the data set:
%
%   library(mixsqp)
%   set.seed(1)
%   n <- 1e5
%   m <- 10
%   L <- simulatemixdata(1e5,10)$L
%   write.table(format(L,digits = 4,scientific = TRUE),
%               "simdata.csv",sep = ",",quote = FALSE,
%               col.names = FALSE,row.names = FALSE)
%
% The solution is, roughly:
%
%   0.492
%   0.000
%   0.000
%   0.000
%   0.006
%   0.161
%   0.333
%   0.009
%   0.000
%   0.000
%

% LOAD DATA
% ---------
% Uncomment this code for other examples:
% load('../analysis/d1.mat');
% L = rand(1000,10000);
% L = csvread('simdata.csv');
L = csvread('simdata-n=1000-m=2000.csv');
[n m] = size(L);

% FIT MODEL
% ---------
x0 = ones(m,1)/m;
fx = @(x) mixobj(L,x,eps);
g  = @(x) projectSimplex(x);
[x f nfevals nproj timings] = ...
    minConf_SPG(fx,x0,g,struct('useSpectral',false,'optTol',1e-8,...
                               'progTol',1e-15,'suffDec',1e-8,...
                               'memory',1,'maxIter',1e4,'verbose',0,...
                               'interp',2,'testOpt',1,'curvilinear',0));

% SAVE RESULTS TO FILE
% --------------------
csvwrite('pg-n=1000-m=2000.csv',[timings f]);