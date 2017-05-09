function [z,dz,ymu,ys,fmu,fs,fpi] = acqThompson(xi,target,gpstruct,optimState,grad)
%ACQTHOMPSON Acquisition function for Thompson sampling (not supported).

if nargin < 5 || isempty(grad); grad = 0; end

n = size(xi,1);
Nhyp = numel(gpstruct.hyp);

if grad == 1 && n > 1
    error('acqNegPI:gradient', ...
        'Gradient of acquisition function is provided only at one test point XI (row vector).');
end

if grad
    [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2,hypw] = gppred(xi,gpstruct);
end
fs = sqrt(fs2);
ys = sqrt(ys2);

post = gpstruct.post;
alpha = post.alpha; L = post.L; sW = post.sW;

Ks = feval(gpstruct.cov{:}, gpstruct.hyp.cov, gpstruct.x, xi);
Kss = feval(gpstruct.cov{:}, gpstruct.hyp.cov, xi);

V  = L'\bsxfun(@times,sW,Ks);
Fs2 = Kss - V'*V;                       % predictive variances
try
    Fs = chol(Fs2);
    z = bsxfun(@plus, fmu(:), Fs'*randn(n,1))';
catch
    z = fmu(:) + fs(:).*randn(n,1);
end
    
    
dz = [];
fpi = [];

end