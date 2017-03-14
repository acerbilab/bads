function [x,fval,exitflag,output] = minimizebnd(fun,x0,LB,UB,len,verbose,varargin)
%MINIMIZEBND Multidimensional nonlinear minimization using conjugate gradients. 
%   X = MINIMIZEBND(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X 
%   and returns a scalar function value F evaluated at X and a vector DF, 
%   the gradient of the function evaluated at X.
%
%   X = MINIMIZEBND(FUN,X0,LB,UB) defines a set of lower and upper bounds 
%   on the design variables, X, so that a solution is found in the range 
%   LB <= X <= UB. Use empty matrices for LB and UB if no bounds exist. Set 
%   LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if X(i) is 
%   unbounded above.
%
%   X = MINIMIZEBND(FUN,X0,LB,UB,LEN) specifies the length of the run LEN.
%   If LEN is positive, it gives the maximum number of line searches, if 
%   negative its absolute gives the maximum allowed number of function 
%   evaluations. The default value is -200*NumberOfVariables. 
%   Optionally, LEN can have a second component, which will indicate the 
%   reduction in function value to be expected in the first line-search
%   (defaults to 1.0).
%
%   X = MINIMIZEBND(FUN,X0,LB,UB,LEN,VERBOSE) specifies whether to display
%   iterations of the algorithm on screen (defaults to 1).
%
%   X = MINIMIZEBND(FUN,X0,LB,UB,LEN,VERBOSE,VARARGIN) specifies additional
%   arguments that are passed to FUN.
%
%   [X,FVAL]= MINIMIZEBND(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = MINIMIZEBND(...) returns an EXITFLAG that describes 
%   the exit condition of MINIMIZEBND. This is kept for compatibility with
%   other MATLAB functions. The default value is 0.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = MINIMIZEBND(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
% Notes:
%
%   MINIMIZEBND is based on FMINSEARCHBND by John D'Errico and MINIMIZE by
%   Carl Edward Rasmussen.
%
%   Variables which are constrained by both a lower and an upper bound will 
%   use a sin transformation. Those constrained by only a lower or an upper 
%   bound will use a quadratic transformation, and unconstrained variables 
%   will be left alone.
%
%   Variables may be fixed by setting their respective bounds equal.
%   In this case, the problem will be reduced in size for MINIMIZE.
%
%   The bounds are inclusive inequalities, which admit the boundary values 
%   themselves, but will not permit ANY function evaluations outside the 
%   bounds. These constraints are strictly followed.
%
%   If your problem has an EXCLUSIVE (strict) constraint which will
%   not admit evaluation at the bound itself, then you must provide
%   a slightly offset bound. An example of this is a function which
%   contains the log of one of its parameters. If you constrain the
%   variable to have a lower bound of zero, then MINIMIZEBND may
%   try to evaluate the function exactly at zero.
%
%   See also FMINSEARCHBND.
%
%
%   Author: Luigi Acerbi
%   E-mail: luigi.acerbi@gmail.com
%   Release date: 07/Oct/2015

% Size checks
xsize = size(x0);
x0 = x0(:);
nvars = numel(x0);

if (nargin < 3) || isempty(LB); LB = -Inf(nvars,1); end
if (nargin < 4) || isempty(UB); UB = Inf(nvars,1); end
LB = LB(:);
UB = UB(:);

if (nvars~=numel(LB)) || (nvars~=numel(UB))
  error ('x0 is incompatible in size with either LB or UB.');
end

% Set default options if necessary
if nargin < 5 || isempty(len); len = -200*nvars; end
if nargin < 6 || isempty(verbose); verbose = 1; end

% stuff into a struct to pass around
params.args = varargin;
params.LB = LB;
params.UB = UB;
params.fun = fun;
params.nvars = nvars;
% note that the number of parameters may actually vary if 
% a user has chosen to fix one or more parameters
params.xsize = xsize;
params.OutputFcn = [];

% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
params.BoundClass = zeros(nvars,1);
k = isfinite(LB) + 2*isfinite(UB);
params.BoundClass = k;
params.BoundClass(k == 3 & LB == UB) = 4;

% Transform starting values into their unconstrained surrogates.
% Check for infeasible starting guesses.
x0u = x0;
k = 1;
for i = 1:nvars
  switch params.BoundClass(i)
    case 1
        % lower bound only
        if x0(i) <= LB(i)
            % infeasible starting value. Use bound.
            x0u(k) = 0;
        else
            x0u(k) = sqrt(x0(i) - LB(i));
        end      
        k=k+1;
      
    case 2
        % Upper bound only
        if x0(i) >= UB(i)
            % infeasible starting value. use bound.
            x0u(k) = 0;
        else
            x0u(k) = sqrt(UB(i) - x0(i));
        end      
        k=k+1;
      
    case 3
        % Lower and upper bounds
        if x0(i) <= LB(i)
            % infeasible starting value
            x0u(k) = -pi/2;
        elseif x0(i) >= UB(i)
            % infeasible starting value
            x0u(k) = pi/2;
        else
            x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
            % shift by 2*pi to avoid problems at zero
            x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
        end      
        k=k+1;
      
    case 0
        % Unconstrained variable. x0u(i) is set.
        x0u(k) = x0(i);      
        k=k+1;
      
    case 4
        % Fixed variable. drop it before minimization sees it.
        % k is not incremented for this variable.
  end
  
end

% If any of the unknowns were fixed, then we need to shorten
% x0u now.
if k <= nvars; x0u(k:nvars) = []; end

% Were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. Quit immediately, setting the
  % appropriate parameters, then return.
  
  % Undo the variable transformations into the original space
  x = xtransform(x0u,params);
  
  % Final reshape
  x = reshape(x,xsize);
  
  % Stuff fval with the final value
  fval = feval(params.fun,x,params.args{:});
  
  % Minimization function was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcount = 1;
  output.algorithm = 'minimize';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to minimize
  return;
end

if all(params.BoundClass == 0) 
    % Pure unconstrained minimization
    [x,fval,iter,funccount] = minimize(x0,fun,len,verbose,params.args{:});
else
    % Constrained minimization with our own intra-objective function.
    [xu,fval,iter,funccount] = minimize(x0u,@intrafun,len,verbose,params);
    
    % Undo the variable transformations into the original space
    x = xtransform(xu,params);
end

exitflag = 0;
output.iterations = iter;
output.funcCount = funccount;
output.algorithm = 'minimize';
if output.funcCount == -len
    output.message = 'Reached maximum number of function evaluations.';
elseif output.iterations == len
    output.message = 'Reached maximum number of iterations (line searches).';
else
    output.message = 'Found local function minimizer.';
end

% Final reshape to make sure the result has the proper shape
x = reshape(x,xsize);

end % mainline end

%--------------------------------------------------------------------------
function [fval,df] = intrafun(x,params)
%INTRAFUN Transform variables, then call original function

% Transform
xtrans = xtransform(x,params);

% Call fun
[fval,df] = feval(params.fun,reshape(xtrans,params.xsize),params.args{:});

% Correct gradient
for k = 1:numel(x)
    switch params.BoundClass(k)
        case 1  % Lower bound only
            df(k) = df(k)*2*x(k);

        case 2  % Upper bound only
            df(k) = -df(k)*2*x(k);

        case 3  % Lower and upper bounds
            df(k) = df(k)*cos(x(k))*(params.UB(k)-params.LB(k))/2;        
    end
end

end % sub function intrafun end

%--------------------------------------------------------------------------
function xtrans = xtransform(x,params)
% Converts unconstrained variables into their original domains

xtrans = zeros(params.xsize);
% k allows some variables to be fixed, thus dropped from the
% optimization.
k=1;
for i = 1:params.nvars
  switch params.BoundClass(i)
    case 1
      % lower bound only
      xtrans(i) = params.LB(i) + x(k).^2;
      
      k=k+1;
    case 2
      % upper bound only
      xtrans(i) = params.UB(i) - x(k).^2;
      
      k=k+1;
    case 3
      % lower and upper bounds
      xtrans(i) = (sin(x(k))+1)/2;
      xtrans(i) = xtrans(i)*(params.UB(i) - params.LB(i)) + params.LB(i);
      % just in case of any floating point problems
      xtrans(i) = max(params.LB(i),min(params.UB(i),xtrans(i)));
      
      k=k+1;
    case 4
      % fixed variable, bounds are equal, set it at either bound
      xtrans(i) = params.LB(i);
    case 0
      % unconstrained variable.
      xtrans(i) = x(k);
      
      k=k+1;
  end
end

end % sub function xtransform end

%--------------------------------------------------------------------------
function [x,fval,iter,funccount] = minimize(x,f,length,verbose,varargin)

% Minimize a differentiable multivariate function using conjugate gradients.
%
% Usage: [X,fval,i] = minimize(X,f,length,verbose,P1,P2,P3, ... )
% 
% X       initial guess
% f       the name or pointer to the function to be minimized. The function
%         f must return two arguments, the value of the function, and its
%         partial derivatives wrt the elements of X. The partial derivative  
%         must have the same type as X.
% length  length of the run; if it is positive, it gives the maximum number of
%         line searches, if negative its absolute gives the maximum allowed
%         number of function evaluations. Optionally, length can have a second
%         component, which will indicate the reduction in function value to be
%         expected in the first line-search (defaults to 1.0).
% verbose verbosity flag; if 1 write iterations, if 0 no output.
% P1, P2, ... parameters are passed to the function f.
%
% X         the returned solution
% fval      function value at solution
% iter      number of iterations (line searches)
% funccount number of function evaluations
%
% The function returns when either its length is up, or if no further progress
% can be made (ie, we are at a (local) minimum, or so close that due to
% numerical problems, we cannot get any closer). NOTE: If the function
% terminates within a few iterations, it could be an indication that the
% function values and derivatives are not consistent (ie, there may be a bug in
% the implementation of your "f" function).
%
% The Polack-Ribiere flavour of conjugate gradients is used to compute search
% directions, and a line search using quadratic and cubic polynomial
% approximations and the Wolfe-Powell stopping criteria is used together with
% the slope ratio method for guessing initial step sizes. Additionally a bunch
% of checks are made to make sure that exploration is taking place and that
% extrapolation will not be unboundedly large.
%
% Copyright (C) 2001 - 2010 by Carl Edward Rasmussen, 2010-01-03.
% Minor interface modifications by Luigi Acerbi 07/Oct/2015.

INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
EXT = 3.0;                  % extrapolate maximum 3 times the current step-size
MAX = 20;                         % max 20 function evaluations per line search
RATIO = 10;                                       % maximum allowed slope ratio
SIG = 0.1; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-
% Powell conditions. SIG is the maximum allowed absolute ratio between
% previous and new slopes (derivatives in the search direction), thus setting
% SIG to low (positive) values forces higher precision in the line-searches.
% RHO is the minimum allowed fraction of the expected (from the slope at the
% initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
% Tuning of SIG (depending on the nature of the function to be optimized) may
% speed up the minimization; it is probably not worth playing much with RHO.

% The code falls naturally into 3 parts, after the initial line search is
% started in the direction of steepest descent. 1) we first enter a while loop
% which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
% have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
% enter the second loop which takes p2, p3 and p4 chooses the subinterval
% containing a (local) minimum, and interpolates it, unil an acceptable point
% is found (Wolfe-Powell conditions). Note, that points are always maintained
% in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
% conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
% was a problem in the previous line-search. Return the best value so far, if
% two consecutive line-searches fail, or whenever we run out of function
% evaluations or line-searches. During extrapolation, the "f" function may fail
% either with an error or returning Nan or Inf, and minimize should handle this
% gracefully.

if nargin < 4; verbose = 0; end

if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
if length>0, S='Linesearch'; else S='Function evaluation'; end 

i = 0;                                            % zero the run length counter
iter = 0;
funccount = 1;
ls_failed = 0;                             % no previous line search has failed
[f0,df0] = feval(f, x, varargin{:});          % get function value and gradient
Z = x;
if verbose > 0; fprintf('%s %6i;  Value %4.6e\r', S, i, f0); end
if exist('fflush','builtin') fflush(stdout); end
fX = f0;
i = i + (length<0);                                            % count epochs?!
s = -df0; d0 = -s'*s;           % initial search direction (steepest) and slope
x3 = red/(1-d0);                                  % initial step is red/(|s|+1)

while i < abs(length)                                      % while not finished
  i = i + (length>0);                                      % count iterations?!
  iter = iter + 1;

  X0 = x; F0 = f0; dF0 = df0;                   % make a copy of current values
  if length>0, M = MAX; else M = min(MAX, -length-i); end

  while 1                             % keep extrapolating as long as necessary
    x2 = 0; f2 = f0; d2 = d0; f3 = f0; df3 = df0;
    success = 0;
    while ~success && M > 0
      try
        M = M - 1; i = i + (length<0); funccount = funccount + 1; % count epochs
        
        [f3,df3] = feval(f, x+x3*s, varargin{:});
        if isnan(f3) || isinf(f3) || any(isnan(df3)+isinf(df3)), error(' '),end
        success = 1;
      catch                                % catch any error which occured in f
        x3 = (x2+x3)/2;                                  % bisect and try again
      end
    end
    if f3 < F0, X0 = x+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    d3 = df3'*s;                                                    % new slope
    if d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0  % are we done extrapolating?
      break
    end
    x1 = x2; f1 = f2; d1 = d2;                        % move point 2 to point 1
    x2 = x3; f2 = f3; d2 = d3;                        % move point 3 to point 2
    A = 6*(f1-f2)+3*(d2+d1)*(x2-x1);                 % make cubic extrapolation
    B = 3*(f2-f1)-(2*d1+d2)*(x2-x1);
    x3 = x1-d1*(x2-x1)^2/(B+sqrt(B*B-A*d1*(x2-x1))); % num. error possible, ok!
    if ~isreal(x3) || isnan(x3) || isinf(x3) || x3 < 0 % num prob | wrong sign?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 > x2*EXT                  % new point beyond extrapolation limit?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 < x2+INT*(x2-x1)         % new point too close to previous point?
      x3 = x2+INT*(x2-x1);
    end
  end                                                       % end extrapolation

  while (abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0  % keep interpolating
    if d3 > 0 || f3 > f0+x3*RHO*d0                         % choose subinterval
      x4 = x3; f4 = f3; d4 = d3;                      % move point 3 to point 4
    else
      x2 = x3; f2 = f3; d2 = d3;                      % move point 3 to point 2
    end
    if f4 > f0           
      x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2));  % quadratic interpolation
    else
      A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);                    % cubic interpolation
      B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
      x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A;        % num. error possible, ok!
    end
    if isnan(x3) || isinf(x3)
      x3 = (x2+x4)/2;               % if we had a numerical problem then bisect
    end
    x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));  % don't accept too close
    [f3,df3] = feval(f, rewrap(Z,x+x3*s), varargin{:});
    if f3 < F0, X0 = x+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    M = M - 1; i = i + (length<0); funccount = funccount + 1;  % count epochs?!
    d3 = df3'*s;                                                    % new slope
  end                                                       % end interpolation

  if abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0          % if line search succeeded
    x = x+x3*s; f0 = f3; fX = [fX' f0]';                     % update variables
    if verbose > 0; fprintf('%s %6i;  Value %4.6e\r', S, i, f0); end
    if exist('fflush','builtin') fflush(stdout); end
    s = (df3'*df3-df0'*df3)/(df0'*df0)*s - df3;   % Polack-Ribiere CG direction
    df0 = df3;                                               % swap derivatives
    d3 = d0; d0 = df0'*s;
    if d0 > 0                                      % new slope must be negative
      s = -df0; d0 = -s'*s;                  % otherwise use steepest direction
    end
    x3 = x3 * min(RATIO, d3/(d0-realmin));          % slope ratio but max RATIO
    ls_failed = 0;                              % this line search did not fail
  else
    x = X0; f0 = F0; df0 = dF0;                     % restore best point so far
    if ls_failed || i > abs(length)         % line search failed twice in a row
      break;                             % or we ran out of time, so we give up
    end
    s = -df0; d0 = -s'*s;                                        % try steepest
    x3 = 1/(1-d0);                     
    ls_failed = 1;                                    % this line search failed
  end
end

if verbose > 0
    fprintf('\n'); if exist('fflush','builtin'); fflush(stdout); end;
end

fval = fX(end);

end

