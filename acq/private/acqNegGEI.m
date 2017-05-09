function [z,dz,ymu,ys,fmu,fs,fpi] = acqNegGEI(xi,target,gpstruct,optimState,grad)
%ACQNEGEI Acquisition function for (negative) generalized expected improvement (unsupported).

if nargin < 5 || isempty(grad); grad = 0; end

n = size(xi,1);

if grad == 1 && n > 1
    error('acqNegGEI:gradient', ...
        'Gradient of acquisition function is provided only at one test point XI (row vector).');
end

if grad
    [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2,hypw] = gppred(xi,gpstruct);
end
fs = sqrt(fs2);
ys = sqrt(ys2);

% Probability of improvement
gammaz = real((target - fmu)./fs);
fpi = 0.5*erfc(-gammaz/sqrt(2));            

% Squared expected improvement
z1 = -(fs.*(gammaz.*fpi + exp(-0.5*(gammaz.^2))/sqrt(2*pi))).^2;

% Expected squared improvement
z2 = -fs.^2.*((gammaz.^2+1).*fpi + gammaz.*exp(-0.5*(gammaz.^2))/sqrt(2*pi));

% gamma(1,1,:) = linspace(0,1,11);
gamma(1,1,:) = 0.5;
z = bsxfun(@times,gamma,z1) + bsxfun(@times,1-gamma,z2);

z = sum(bsxfun(@times,hypw(~isnan(hypw)),z(~isnan(hypw),:,:)),1);

[zi,idx] = min(z,[],3);
[~,idx2] = min(zi,[],2);
idx = idx(idx2);
z = z(:,:,idx);

%try
%    z = sum(bsxfun(@times,hypw(~isnan(hypw)),z(~isnan(hypw),:)),1);
%catch
%    z = Inf(1,n);
%    dz = NaN(n,size(xi,2));
%    return;
%end

if grad    
    % Gradient of probability of improvement
    dfs = 0.5*dfs2./fs;
    dgammaz = -(dfmu.*fs + (target - fmu).*dfs)./fs2;
    dfpi = -0.5*dgammaz/sqrt(2)*(-2*exp(-gammaz.^2/2)/sqrt(pi));
    
    % Gradient of expected improvement
    %dz = -(dfs.*gammaz.*fpi + dgammaz.*fs.*fpi + dfpi.*gammaz.*fs) ...
    %    + (fs.*gammaz.*dgammaz - dfs).*exp(-0.5*gammaz.^2)/sqrt(2*pi);
    %dz = sum(bsxfun(@times,hypw,dz(~isnan(hypw),:)),1);    
else
    dz = NaN(n,size(xi,2));     % Gradient not estimated
end

end