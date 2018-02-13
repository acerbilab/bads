function [hyp,exitflag] = gpHyperOptimize(hyp0,gpstruct,optimtoolbox,options,nudge,removeafter)
%GPHYPEROPTIMIZE Optimize Gaussian Process hyperparameters.
      
if nargin < 3; optimtoolbox = 1; end
if nargin < 4; options = []; end
if nargin < 5 || isempty(nudge); nudge = [0 0]; end
if numel(nudge) == 1; nudge(2) = 0.5*nudge(1); end
if nargin < 6 || isempty(removeafter); removeafter = Inf; end

% Perform initial argument checks/transformations
gpstruct = gpset(gpstruct);

nvars = size(gpstruct.x, 2);

% Find fixed parameters (delta or clamped priors), they are not optimized
fixed = []; offset = 0;
for field = {'cov', 'lik', 'mean'}
    num_entries = numel(gpstruct.prior.(field{:}));
    fixed = [fixed,zeros(1,num_entries)];
    for j = 1:num_entries
        if isequal(gpstruct.prior.(field{:}){j}{1},@delta_prior) || ...
                strcmp(gpstruct.prior.(field{:}){j}{1}, 'delta_prior') || ...
                isequal(gpstruct.prior.(field{:}){j}{1},@priorClamped) || ...
                strcmp(gpstruct.prior.(field{:}){j}{1}, 'priorClamped') || ...
                isequal(gpstruct.prior.(field{:}){j}{1},@priorDelta) || ...
                strcmp(gpstruct.prior.(field{:}){j}{1}, 'priorDelta')
            fixed(offset+j) = 1;
        end
    end
    offset = offset + num_entries;      
end
fixed = logical(fixed);

% Get bounds
bounds = unwrap2vec(gpstruct.bounds);
LB = bounds(1:2:end-1); UB = bounds(2:2:end);
lb = LB(~fixed); ub = UB(~fixed);

[func,hessfunc,optfunc,OptOptions] = getOptimizer(lb,ub,fixed,gpstruct,optimtoolbox,options);

fval = zeros(1,numel(hyp0));
theta = zeros(sum(~fixed),numel(hyp0));

Nruns = numel(hyp0);
success_vec = zeros(1,Nruns);

for iRun = 1:Nruns
    if iscell(hyp0); hyp = hyp0{iRun}; else hyp = hyp0(iRun); end
    
    % Unwrap starting point, force it to be inside bounds
    theta0 = unwrap2vec(hyp);
    theta0 = min(max(theta0(~fixed),lb),ub);
        
    % Evaluate starting point
    try
        f0 = func(theta0);
    catch
        f0 = Inf;
    end
    
    theta(:,iRun) = theta0;
    fval(iRun) = f0;
   
    nTry = 10;     % Try up to ten times    
    noiseNudge = 0;
              
    for iTry = 1:nTry
        % iTry
        
        % Require a minimum number of points to do the fit
        if size(gpstruct.y,1) < nvars; break; end
        
        try
            tic
            [theta(:,iRun),fval(iRun)] = optfunc(theta0,OptOptions);
            t1 = toc;

            if 0
                gpstruct.fixed = fixed;
                LB = lb(~fixed); LB(isinf(LB)) = -10;
                UB = ub(~fixed); UB(isinf(UB)) = 10;
                iinit.x0 = [LB, theta0(:), UB];
                iinit.l = 2*ones(numel(theta0),1);
                iinit.L = 3*ones(numel(theta0),1);
                smax = 5*numel(theta0)+10;
                stop = 3*numel(theta0);

                % optfunc2 = @(x0_,opt_) mymcs('mcs_gp_optimizer',gpstruct,lb(~fixed),ub(~fixed),1,smax,500,stop,iinit);
                OptOptions = optimset(options,'Display','off','GradObj','on');
                if strcmp(OptOptions.Hessian,'user-supplied')
                    OptOptions = optimset(OptOptions,'HessFcn',hessfunc);
                end
                optfunc2 = @(x0_,opt_) fmincon(func,x0_,[],[],[],[],lb(~fixed),ub(~fixed),[],opt_);

                tic
                [xbest2,fval2] = optfunc2(theta0,OptOptions);
                t2 = toc;
                [t1 t2]
                [fval(iRun) fval2]
            end
            
            success_vec(iRun) = 1;

            break;
        catch optexc
            
            
            % Try removing closest pair of points
            if iTry > removeafter % if iTry > 2
                
                %optexc
                
                % gpstruct
                idx = false(size(gpstruct.y));
                if 1
                    x1 = []; x2 = [];

                    x1(:,:,1) = gpstruct.x;
                    x2(1,:,:) = gpstruct.x';
                    dist = squeeze(sum(bsxfun(@minus, x1, x2).^2, 2));
                    dist = dist + tril(Inf(size(dist)));
                    [mins,row] = min(dist,[],1);
                    [d2,col] = min(mins,[],2);
                    row = row(col);
                    % sqrt(d2)


                    % gpstruct.hyp.cov

                    if gpstruct.y(row) > gpstruct.y(col)
                        idx(row) = true;
                    else
                        idx(col) = true;
                    end                    
                    % [~,idx] = max(gpstruct.y);
                end
                
                idx = idx | gpstruct.y > prctile1(gpstruct.y,95);
                
                gpstruct.x(idx,:) = [];
                gpstruct.y(idx) = [];
                                
                % Recompute optimization function
                [func,hessfunc,optfunc,OptOptions] = getOptimizer(lb,ub,fixed,gpstruct,optimtoolbox,options);
            end
            
            % Retry with random sample from prior        
            theta0 = unwrap2vec(gppriorrnd(gpstruct.prior,gpstruct.hyp(1)));
            
            % Mean parameter
            %mu = median(gpstruct.y);
            %sigma = std(gpstruct.y);
            %theta0(end) = sigma*randn() + mu;
            
            thetaold = unwrap2vec(hyp);
            theta0 = 0.5*(theta0 + thetaold);
            
            % Try increase starting point of noise
            noiseNudge = noiseNudge + nudge(1);
            
            % Increase lower bound on noise
            lb(end-1) = lb(end-1) + nudge(2);
                        
            theta0(end-1) = theta0(end-1)+noiseNudge;
            theta0 = min(max(theta0(~fixed),lb),ub);
            
            [func,hessfunc,optfunc,OptOptions] = getOptimizer(lb,ub,fixed,gpstruct,optimtoolbox,options);
            
            % theta0'
            
            %theta0'
            %warning('Need to fix this!')
            
        end
    end
    
    if ~success_vec(iRun)
        if isfield(options,'gpWarnings') && ~isempty(options.gpWarnings) ...
                && options.gpWarnings
            warning('bads:gpHyperOptFail', ['Failed optimization of hyper-parameters (' num2str(nTry) ' attempts). GP approximation might be unreliable.']);
        end
        
        if 0
            x1(:,:,1) = gpstruct.x;
            x2(1,:,:) = gpstruct.x';
            dist = sqrt(squeeze(sum(bsxfun(@minus, x1, x2).^2, 2)));
            save('error.mat');
            error('a')
        end
        
        
    end
end

[~,index] = min(fval);
newtheta = unwrap2vec(gpstruct.hyp(1));
newtheta(~fixed) = theta(:,index);
hyp = rewrap(gpstruct.hyp(1), newtheta);

% Set EXITFLAG
if all(success_vec)     % All runs succeded
    exitflag = 1;    
elseif any(success_vec) % At least one run succeded
    exitflag = 0;
else
    exitflag = -1;      % All runs failed
end

% exp(hyp.cov(1:size(gpstruct.x,2)))'
%hyp.cov(:)'
%hyp.lik(:)'
%hyp.mean

%figure(100);
%landscapeplot(func,theta(:,index)',lb(~fixed),ub(~fixed),0.5,[]);
%drawnow;

%hlZ = hessian_wrapper(theta(:,index),gpstruct,fixed);
%[v,e] = eig(inv(hlZ));
%[sqrt(e), v]

end

%--------------------------------------------------------------------------
%GETOPTIMIZED Choose optimizer -- see if Optimization Toolbox is available
function [func,hessfunc,optfunc,OptOptions] = getOptimizer(lb,ub,fixed,gpstruct,optimtoolbox,options)

% Constrained or unconstrained optimization?
if any(lb > -Inf) || any(ub < Inf); ConOpt = 1; else ConOpt = 0; end

% Check if MATLAB's Optimization Toolbox is available
if isempty(optimtoolbox)
    if exist('fmincon.m','file') && exist('fminunc.m','file') && exist('optimoptions.m','file')
        optimtoolbox = 1;
    else
        optimtoolbox = 0;
    end    
end

% Optimization functions
func = @(theta_) gp_optimizer(theta_, fixed, gpstruct);
hessfunc = @(theta_,lambda_) hessian_wrapper(theta_,gpstruct,fixed);

% Number of variables for optimization
numberOfVariables = sum(~fixed);

if optimtoolbox
    if ConOpt
        % Constrained problem -- use FMINCON
        OptOptions = optimset(options,'Display','off','GradObj','on','DerivativeCheck','off','Algorithm','interior-point');
        if strcmp(OptOptions.Hessian,'user-supplied')
            OptOptions = optimset(OptOptions,'HessFcn',hessfunc);
        end
        optfunc = @(x0_,opt_) fmincon(func,x0_,[],[],[],[],lb,ub,[],opt_);
    else
        % Unconstrained problem -- use FMINUNC 
        OptOptions = optimset(options,'Display','off','GradObj','on');
        if strcmp(OptOptions.Hessian,'user-supplied')
            OptOptions = optimset(OptOptions,'HessFcn',hessfunc);
        end
        optfunc = @(x0_,opt_) fminunc(func,x0_,opt_);
    end
else
    % Optimization Toolbox is not available, use MINIMIZEBND
    if isfield(options,'MaxFunEval') && ~isempty(options.MaxFunEval)
        len = -abs(options.MaxFunEval);
    elseif isfield(options,'MaxFunEvals') && ~isempty(options.MaxFunEvals)
        len = -abs(options.MaxFunEvals);
    elseif isfield(options,'MaxIter') && ~isempty(options.MaxIter)
        len = abs(options.MaxIter);
    else
        len = -100*numberOfVariables;
    end
    if isfield(options,'Display') && ~isempty(options.Display)
        verbose = ~isempty(strfind(lower(options.Display),'iter'));
    else
        verbose = 0;
    end
    OptOptions = [];
    optfunc = @(x0_,opt_) minimizebnd(func,x0_,lb,ub,len,verbose);    
end

end

%--------------------------------------------------------------------------
function [nlZ, dnlZ, HnlZ] = gp_optimizer(theta, fixed, gpstruct)
%GP_OPTIMIZER Wrapper function for GP optimization

newtheta = unwrap2vec(gpstruct.hyp(1));
newtheta(~fixed) = theta;
thetastruct = rewrap(gpstruct.hyp(1), newtheta);

f = @() (feval(gpstruct.inf{:}, thetastruct, gpstruct.mean, gpstruct.cov, gpstruct.lik, gpstruct.x, gpstruct.y));

if nargout <= 1
    [~, nlZ] = f();
    return;

elseif nargout == 2
    [~, nlZ, dnlZ] = f();

elseif nargout > 2
    [~, nlZ, dnlZ, ~, ~, HnlZ] = f();

    HnlZ = HnlZ.value;
    HnlZ(fixed,:) = [];
    HnlZ(:,fixed) = [];
end

dnlZ = unwrap2vec(dnlZ);
dnlZ(fixed) = [];
    
end

%--------------------------------------------------------------------------
function [hlZ] = hessian_wrapper(theta,gpstruct,fixed)
%HESSIAN_WRAPPER Compute Hessian for a given hyper-parameter vector THETA

[~,~,hlZ] = gp_optimizer(theta, fixed, gpstruct);        

end

%--------------------------------------------------------------------------
function [x,fval,i] = gpminimize(x,f,length,verbose,varargin)

% Minimize a differentiable multivariate function using conjugate gradients.
%
% Usage: [X,fval,i] = minimize(X,f,length,verbose,P1,P2,P3, ... )
% 
% X       initial guess; may be of any type, including struct and cell array
% f       the name or pointer to the function to be minimized. The function
%         f must return two arguments, the value of the function, and it's
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
% X       the returned solution
% fval    function value at solution
% i       number of iterations (line searches or function evaluations, 
%         depending on the sign of "length") used at termination.
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
% See also: checkgrad 
%
% Copyright (C) 2001 - 2010 by Carl Edward Rasmussen, 2010-01-03

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

  X0 = x; F0 = f0; dF0 = df0;                   % make a copy of current values
  if length>0, M = MAX; else M = min(MAX, -length-i); end

  while 1                             % keep extrapolating as long as necessary
    x2 = 0; f2 = f0; d2 = d0; f3 = f0; df3 = df0;
    success = 0;
    while ~success && M > 0
      try
        M = M - 1; i = i + (length<0);                         % count epochs?!
        
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
    M = M - 1; i = i + (length<0);                             % count epochs?!
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





