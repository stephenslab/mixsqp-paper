function [x, obj, funEvals, projects, timings] = ...
  minConf_SPG(funObj, x, funProj, options)
% function [x,f] = minConF_SPG(funObj,x,funProj,options)
%
% Function for using Spectral Projected Gradient to solve problems of the form
%   min funObj(x) s.t. x in C
%
%   @funObj(x): function to minimize (returns gradient as second argument)
%   @funProj(x): function that returns projection of x onto C
%
%   options:
%       verbose: level of verbosity (0: no output, 1: final, 2: iter (default), 3:
%       debug)
%       optTol: tolerance used to check for optimality (default: 1e-5)
%       progTol: tolerance used to check for lack of progress (default: 1e-9)
%       maxIter: maximum number of calls to funObj (default: 500)
%       numDiff: compute derivatives numerically (0: use user-supplied
%       derivatives (default), 1: use finite differences, 2: use complex
%       differentials)
%       suffDec: sufficient decrease parameter in Armijo condition (default
%       : 1e-4)
%       interp: type of interpolation (0: step-size halving, 1: quadratic,
%       2: cubic)
%       memory: number of steps to look back in non-monotone Armijo
%       condition
%       useSpectral: use spectral scaling of gradient direction (default:
%       1)
%       curvilinear: backtrack along projection Arc (default: 0)
%       testOpt: test optimality condition (default: 1)
%       feasibleInit: if 1, then the initial point is assumed to be
%       feasible
%       bbType: type of Barzilai Borwein step (default: 1)
%
%   Notes: 
%       - if the projection is expensive to compute, you can reduce the
%           number of projections by setting testOpt to 0


nVars = length(x);

% Set Parameters
if nargin < 4
    options = [];
end
[verbose,numDiff,optTol,progTol,maxIter,suffDec,interp,memory,useSpectral,curvilinear,feasibleInit,testOpt,bbType] = ...
    myProcessOptions(...
    options,'verbose',2,'numDiff',0,'optTol',1e-5,'progTol',1e-9,'maxIter',500,'suffDec',1e-4,...
    'interp',2,'memory',10,'useSpectral',1,'curvilinear',0,'feasibleInit',0,...
    'testOpt',1,'bbType',1);

% Output Log
if verbose >= 2
    if testOpt
        fprintf('%10s %10s %10s %15s %15s %14s\n','Iteration','FunEvals','Projections','Step Length','Function Val','Opt Cond');
    else
        fprintf('%10s %10s %10s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val');
    end
end

% Make objective function (if using numerical derivatives)
funEvalMultiplier = 1;
if numDiff
    if numDiff == 2
        useComplex = 1;
    else
        useComplex = 0;
    end
    funObj = @(x)autoGrad(x,useComplex,funObj);
    funEvalMultiplier = nVars+1-useComplex;
end

% Evaluate Initial Point
if ~feasibleInit
    x = funProj(x);
end
[f,g] = funObj(x);
projects = 1;
funEvals = 1;
obj     = zeros(maxIter,1);
timings = zeros(maxIter,1);

% Optionally check optimality
if testOpt
    projects = projects+1;
    if max(abs(funProj(x-g)-x)) < optTol
        if verbose >= 1
        fprintf('First-Order Optimality Conditions Below optTol at Initial Point\n');
        end
        return;
    end
end

i = 1;
while funEvals <= maxIter
    tic;

    % Compute Step Direction
    if i == 1 || ~useSpectral
        alpha = 1;
    else
        y = g-g_old;
        s = x-x_old;
        if bbType == 1
            alpha = (s'*s)/(s'*y);
        else
            alpha = (s'*y)/(y'*y);
        end
        if alpha <= 1e-10 || alpha > 1e10
            alpha = 1;
        end
    end
    d = -alpha*g;
    f_old = f;
    x_old = x;
    g_old = g;
    
    % Compute Projected Step
    if ~curvilinear
      d0 = d;
      d  = funProj(x + d0) - x;
      projects = projects+1;
    end

    % Check that Progress can be made along the direction
    gtd = g'*d;

    % Select Initial Guess to step length
    if i == 1
        t = min(1,1/sum(abs(g)));
    else
        t = 1;
    end

    % Compute reference function for non-monotone condition

    if memory == 1
        funRef = f;
    else
        if i == 1
            old_fvals = repmat(-inf,[memory 1]);
        end

        if i <= memory
            old_fvals(i) = f;
        else
            old_fvals = [old_fvals(2:end);f];
        end
        funRef = max(old_fvals);
    end

    % Evaluate the Objective and Gradient at the Initial Step
    if curvilinear
        x_new = funProj(x + t*d);
        projects = projects+1;
    else
        x_new = x + t*d;
    end
    [f_new,g_new] = funObj(x_new);
    funEvals = funEvals+1;

    % Backtracking Line Search
    lineSearchIters = 1;
    while f_new > funRef + suffDec*g'*(x_new-x) || ~isLegal(f_new)
        temp = t;
        if interp == 0 || ~isLegal(f_new)
            if verbose == 3
                fprintf('Halving Step Size\n');
            end
            t = t/2;
        elseif interp == 2 && isLegal(g_new)
            if verbose == 3
                fprintf('Cubic Backtracking\n');
            end
            t = polyinterp([0 f gtd; t f_new g_new'*d]);
        elseif lineSearchIters < 2 || ~isLegal(f_prev)
            if verbose == 3
                fprintf('Quadratic Backtracking\n');
            end
            t = polyinterp([0 f gtd; t f_new sqrt(-1)]);
        else
            if verbose == 3
                fprintf('Cubic Backtracking on Function Values\n');
            end
            t = polyinterp([0 f gtd; t f_new sqrt(-1);t_prev f_prev sqrt(-1)]);
        end

        % Adjust if change is too small
        if t < temp*1e-3
            if verbose == 3
                fprintf('Interpolated value too small, Adjusting\n');
            end
            t = temp*1e-3;
        elseif t > temp*0.6
            if verbose == 3
                fprintf('Interpolated value too large, Adjusting\n');
            end
            t = temp*0.6;
        end

        % Check whether step has become too small
        if max(abs(t*d)) < progTol || t == 0
            if verbose == 3
                fprintf('Line Search failed\n');
            end
            t = 0;
            f_new = f;
            g_new = g;
            break;
        end

        % Evaluate New Point
        f_prev = f_new;
        t_prev = temp;
        if curvilinear
            x_new = funProj(x + t*d);
            projects = projects+1;
        else
            x_new = x + t*d;
        end
        [f_new,g_new] = funObj(x_new);
        funEvals = funEvals+1;
        lineSearchIters = lineSearchIters+1;

    end

    % Take Step
    x = x_new;
    f = f_new;
    g = g_new;

    if testOpt
        optCond = max(abs(funProj(x-g)-x));
        projects = projects+1;
    end

    % Output Log
    if verbose >= 2
        if testOpt
            fprintf('%10d %10d %10d %15.5e %15.9e %14.8e\n',i,funEvals*funEvalMultiplier,projects,t,f,optCond);
        else
            fprintf('%10d %10d %10d %15.5e %15.5e\n',i,funEvals*funEvalMultiplier,projects,t,f);
        end
    end

    % Check optimality
    if testOpt
        if optCond < optTol
            if verbose >= 1
            fprintf('First-Order Optimality Conditions Below optTol\n');
            end
            break;
        end
    end

    if funEvals*funEvalMultiplier > maxIter
        if verbose >= 1
            fprintf('Function Evaluations exceeds maxIter\n');
        end
        break;
    end

    obj(i)     = f;
    timings(i) = toc;
    i = i + 1;
end
obj     = obj(1:(i-1));
timings = timings(1:(i-1));
end

% ----------------------------------------------------------------------------
function [minPos,fmin] = polyinterp(points,doPlot,xminBound,xmaxBound)
% function [minPos] = polyinterp(points,doPlot,xminBound,xmaxBound)
%
%   Minimum of interpolating polynomial based on function and derivative
%   values
%
%   In can also be used for extrapolation if {xmin,xmax} are outside
%   the domain of the points.
%
%   Input:
%       points(pointNum,[x f g])
%       doPlot: set to 1 to plot, default: 0
%       xmin: min value that brackets minimum (default: min of points)
%       xmax: max value that brackets maximum (default: max of points)
%
%   set f or g to sqrt(-1) if they are not known
%   the order of the polynomial is the number of known f and g values minus 1

if nargin < 2
    doPlot = 0;
end

nPoints = size(points,1);
order = sum(sum((imag(points(:,2:3))==0)))-1;

% Code for most common case:
%   - cubic interpolation of 2 points
%       w/ function and derivative values for both
%   - no xminBound/xmaxBound

if nPoints == 2 && order ==3 && nargin <= 2 && doPlot == 0
    % Solution in this case (where x2 is the farthest point):
    %    d1 = g1 + g2 - 3*(f1-f2)/(x1-x2);
    %    d2 = sqrt(d1^2 - g1*g2);
    %    minPos = x2 - (x2 - x1)*((g2 + d2 - d1)/(g2 - g1 + 2*d2));
    %    t_new = min(max(minPos,x1),x2);
    [minVal minPos] = min(points(:,1));
    notMinPos = -minPos+3;
    d1 = points(minPos,3) + points(notMinPos,3) - 3*(points(minPos,2)-points(notMinPos,2))/(points(minPos,1)-points(notMinPos,1));
    d2 = sqrt(d1^2 - points(minPos,3)*points(notMinPos,3));
    if isreal(d2)
        t = points(notMinPos,1) - (points(notMinPos,1) - points(minPos,1))*((points(notMinPos,3) + d2 - d1)/(points(notMinPos,3) - points(minPos,3) + 2*d2));
        minPos = min(max(t,points(minPos,1)),points(notMinPos,1));
    else
        minPos = mean(points(:,1));
    end
    return;
end

xmin = min(points(:,1));
xmax = max(points(:,1));

% Compute Bounds of Interpolation Area
if nargin < 3
    xminBound = xmin;
end
if nargin < 4
    xmaxBound = xmax;
end

% Constraints Based on available Function Values
A = zeros(0,order+1);
b = zeros(0,1);
for i = 1:nPoints
    if imag(points(i,2))==0
        constraint = zeros(1,order+1);
        for j = order:-1:0
            constraint(order-j+1) = points(i,1)^j;
        end
        A = [A;constraint];
        b = [b;points(i,2)];
    end
end

% Constraints based on available Derivatives
for i = 1:nPoints
    if isreal(points(i,3))
        constraint = zeros(1,order+1);
        for j = 1:order
            constraint(j) = (order-j+1)*points(i,1)^(order-j);
        end
        A = [A;constraint];
        b = [b;points(i,3)];
    end
end

% Find interpolating polynomial
params = A\b;

% Compute Critical Points
dParams = zeros(order,1);
for i = 1:length(params)-1
    dParams(i) = params(i)*(order-i+1);
end

if any(isinf(dParams))
    cp = [xminBound;xmaxBound;points(:,1)].';
else
    cp = [xminBound;xmaxBound;points(:,1);roots(dParams)].';
end

% Test Critical Points
fmin = inf;
minPos = (xminBound+xmaxBound)/2; % Default to Bisection if no critical points valid
for xCP = cp
    if imag(xCP)==0 && xCP >= xminBound && xCP <= xmaxBound
        fCP = polyval(params,xCP);
        if imag(fCP)==0 && fCP < fmin
            minPos = real(xCP);
            fmin = real(fCP);
        end
    end
end
% Plot Situation
if doPlot
    figure(1); clf; hold on;

    % Plot Points
    plot(points(:,1),points(:,2),'b*');

    % Plot Derivatives
    for i = 1:nPoints
        if isreal(points(i,3))
            m = points(i,3);
            b = points(i,2) - m*points(i,1);
            plot([points(i,1)-.05 points(i,1)+.05],...
                [(points(i,1)-.05)*m+b (points(i,1)+.05)*m+b],'c.-');
        end
    end

    % Plot Function
    x = min(xmin,xminBound)-.1:(max(xmax,xmaxBound)+.1-min(xmin,xminBound)-.1)/100:max(xmax,xmaxBound)+.1;
    size(x)
    for i = 1:length(x)
        f(i) = polyval(params,x(i));
    end
    plot(x,f,'y');
    axis([x(1)-.1 x(end)+.1 min(f)-.1 max(f)+.1]);

    % Plot Minimum
    plot(minPos,fmin,'g+');
    if doPlot == 1
        pause(1);
    end
end
end

% ----------------------------------------------------------------------------
function [legal] = isLegal(v)
  legal = sum(any(imag(v(:))))==0 & sum(isnan(v(:)))==0 & sum(isinf(v(:)))==0;
end

% ----------------------------------------------------------------------------
function [varargout] = myProcessOptions(options,varargin)
% Similar to processOptions, but case insensitive and
%   using a struct instead of a variable length list

options = toUpper(options);

for i = 1:2:length(varargin)
    if isfield(options,upper(varargin{i}))
        v = getfield(options,upper(varargin{i}));
        if isempty(v)
            varargout{(i+1)/2}=varargin{i+1};
        else
            varargout{(i+1)/2}=v;
        end
    else
        varargout{(i+1)/2}=varargin{i+1};
    end
end
end

% ----------------------------------------------------------------------------
function [o] = toUpper(o)
if ~isempty(o)
    fn = fieldnames(o);
    for i = 1:length(fn)
        o = setfield(o,upper(fn{i}),getfield(o,fn{i}));
    end
end
end