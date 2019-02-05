% A short script to test the projected gradient method on a very small
% maximum-likelihood estimation problem with a 4 x 3 matrix. The solution
% should be [1/3 2/3 0].
%
% CREATE DATA SET
% ---------------
e = 1e-8;

% R code:
%
%   L <- rbind(c(1,e,e),
%              c(e,1,e),
%              c(e,1,e),
%              c(1,1,1))
%
L = [ 1 e e
      e 1 e
      e 1 e
      1 1 1 ];

% FIT MODEL
% ---------
x0 = [1 1 1]'/3;
f  = @(x) mixobj(L,x,eps);
g  = @(x) projectSimplex(x);
x  = minConf_SPG(f,x0,g,struct('useSpectral',false,'optTol',1e-8,...
                               'progTol',1e-15,'suffDec',1e-8,...
                               'memory',1));
