function [z,dz,ymu,ys,fmu,fs,fpi] = acqLCB(sqrtbetat,xi,target,gpstruct,optimState,grad_flag)
%ACQLCB Acquisition function for lowest confidence bound (LCB).

if nargin < 6 || isempty(grad_flag); grad_flag = false; end

n = size(xi,1);
nvars = size(xi,2);
t = optimState.funccount + 1;

if isempty(sqrtbetat)
    % By default use beta_t from Theorem 1 in Srinivas et al. (2010)
    delta = 0.1;
    nu = 0.2;      % Empirical correction
    sqrtbetat = sqrt(nu*2*log(nvars*t^2*pi^2/(6*delta)));
    % sqrtbetat = sqrt(2*log(t^(nvars/2+2)*pi^2/(3*delta)));
elseif ischar(sqrtbetat) || isa(sqrtbetat,'function_handle')
    sqrtbetat = feval(sqrtbetat,t,nvars);
elseif ~isnumeric(sqrtbetat) || numel(sqrtbetat) > 1
    error('acqLCB:boundparam', ...
        'The SQRTBETAT parameter of the acquisition function needs to be a scalar or a function handle/name to an annealing schedule.');
end

if grad_flag
    [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2,hypw] = gppred(xi,gpstruct);
end
fs = sqrt(fs2);
ys = sqrt(ys2);

% Lower confidence bound
z = fmu - sqrtbetat.*fs;

try
    z = sum(bsxfun(@times,hypw(~isnan(hypw)), z(~isnan(hypw),:) ),1);  
catch
    z = Inf(1,n);
    dz = NaN(n,size(xi,2));
    return;
end
    
if grad_flag
    dz = dfmu - 0.5*sqrtbetat.*dfs2./fs;
    dz = sum(bsxfun(@times,hypw(~isnan(hypw)),dz(~isnan(hypw),:),1));    
else
    dz = NaN(n,size(xi,2));     % Gradient not estimated
end

end