function [z,dz,ymu,ys,fmu,fs,fpi] = acqNegPI(xi,target,gpstruct,optimState,grad_flag)
%ACQNEGPI Acquisition function for (negative) probability of improvement.

if nargin < 5 || isempty(grad_flag); grad_flag = false; end

n = size(xi,1);
Nhyp = numel(gpstruct.hyp);

if grad_flag && n > 1
    error('acqNegPI:gradient', ...
        'Gradient of acquisition function is provided only at one test point XI (row vector).');
end

if grad_flag
    [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2,hypw] = gppred(xi,gpstruct);
end
fs = sqrt(fs2);
ys = sqrt(ys2);

% Negative probability of improvement
gammaz = (target - fmu)./fs;
z = -0.5*erfc(-gammaz/sqrt(2));            

try
    z = sum(bsxfun(@times,hypw(~isnan(hypw)),z(~isnan(hypw),:)),1);
catch
    z = Inf(1,n);
    dz = NaN(n,size(xi,2));
    return;
end

if grad_flag
    % Gradient of probability of improvement
    dfs = 0.5*dfs2./fs;
    dgammaz = -(dfmu.*fs + (target - fmu).*dfs)./fs2;
    dz = 0.5*dgammaz/sqrt(2)*(-2*exp(-gammaz.^2/2)/sqrt(pi));
    
    dz = sum(bsxfun(@times,hypw,dz(~isnan(hypw),:)),1);    
else
    dz = NaN(n,size(xi,2));     % Gradient not estimated
end

end