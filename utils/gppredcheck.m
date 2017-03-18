function [h,p] = gppredcheck(gpstats,alpha)
%GPPREDCHECK Check calibration of Gaussian process prediction.

n = gpstats.last;
idx = 1:n;
p = NaN;

% No predictions available, test failed
if n == 0; h = 1; return; end

zscores = (gpstats.fval(idx) - gpstats.ymu(:,idx))./(gpstats.ys(:,idx));

% A NaN means that some prediction went wrong
if any(isnan(zscores)); h = 1; return ; end
    
zscores(isnan(zscores)) = [];
n = numel(zscores);

if n < 3
    mychi2inv = @(x,v) 2*gammaincinv(x,v/2);
    plo = mychi2inv(alpha/2,n);
    phi = mychi2inv(1-alpha/2,n);                
    t = sum(zscores.^2);
    if t < plo || t > phi || any(isnan([plo,phi]))
        h = 1; 
    else
        h = 0;
    end                
else
    [h,p] = swtest(zscores, alpha);
end