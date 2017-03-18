function [z,dz,ymu,ys,fmu,fs,fpi] = acqNegEQI(beta,xi,target,gpstruct,optimState,grad)
%ACQNEGEQI Acquisition function for (negative) expected quantile improvement.

if isempty(beta); beta = 0.841344746068543; end
if nargin < 6 || isempty(grad); grad = 0; end

n = size(xi,1);
c = -sqrt(2)*erfcinv(2*beta);

if grad == 1 && n > 1
    error('acqNegEQI:gradient', ...
        'Gradient of acquisition function is provided only at one test point XI (row vector).');
end


if grad
    [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2,hypw] = gppred(xi,gpstruct);
end
fs = sqrt(fs2);
ys = sqrt(ys2);

% Collect observation noise
Nhyp = numel(gpstruct.hyp);
noise = zeros(Nhyp,n);
for i = 1:numel(gpstruct.hyp)
    noise(i,:) = gpstruct.hyp.lik;
end

target = optimState.ftargetmu + c*optimState.ftargets;
qfmu = fmu + c*sqrt(fs2.*noise.^2./(fs2 + noise.^2));
qfs = sqrt(fs2.^2./(fs2 + noise.^2));

% Probability of improvement
qgammaz = real((target - qfmu)./qfs);
qfpi = 0.5*erfc(-qgammaz/sqrt(2));            

% Expected quantile improvement
z = -qfs.*(qgammaz.*qfpi + exp(-0.5*(qgammaz.^2))/sqrt(2*pi));

try
    z = sum(bsxfun(@times,hypw(~isnan(hypw)),z(~isnan(hypw),:)),1);
catch
    z = Inf(1,n);
    dz = NaN(n,size(xi,2));
    return;
end

if grad    
    error('Gradient not supported yet.');
    
    % Gradient of probability of improvement
    dfs = 0.5*dfs2./fs;
    dgammaz = -(dfmu.*fs + (target - fmu).*dfs)./fs2;
    dfpi = -0.5*dgammaz/sqrt(2)*(-2*exp(-gammaz.^2/2)/sqrt(pi));
    
    % Gradient of expected improvement
    dz = -(dfs.*gammaz.*fpi + dgammaz.*fs.*fpi + dfpi.*gammaz.*fs) ...
        + (fs.*gammaz.*dgammaz - dfs).*exp(-0.5*gammaz.^2)/sqrt(2*pi);
    dz = sum(bsxfun(@times,hypw,dz(~isnan(hypw),:)),1);    
else
    dz = NaN(n,size(xi,2));     % Gradient not estimated
end

end