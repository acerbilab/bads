function [z,dz,ymu,ys,fmu,fs,fpi] = acqRnd(xi,target,gpstruct,optimState,grad)
%ACQRND Random acquisition function (not for real usage, for testing only).

if nargin < 5 || isempty(grad); grad = 0; end

n = size(xi,1);
Nhyp = numel(gpstruct.hyp);

if grad == 1 && n > 1
    error('acqRnd:gradient', ...
        'Gradient of acquisition function is provided only at one test point XI (row vector).');
end

ymu = zeros(Nhyp,n);        % Observed value SD
ys = zeros(Nhyp,n);         % Observed value SD
fmu = zeros(Nhyp,n);        % Function mean
fs = zeros(Nhyp,n);         % Function SD

z = randn(1,n);
dz = randn(1,n);

end