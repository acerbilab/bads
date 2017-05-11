function y = prctile1(x,p)
%PRCTILE1 Percentiles of a vector sample.
%   Y = PRCTILE1(X,P) returns percentiles of the values in X.  P is a single
%   percent value, X is a vector, and Y is the P-th percentile.
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,75); % upper IQR
%
%   See also MEDIAN.

if ~isvector(x)
    error('PRCTILE1 requires X to be a vector.');
end

if ~isscalar(p) || p < 0 || p > 100
    error('PRCTILE1 requires P to be a scalar between 0 and 100.');    
end

% Prepare input and remove NaNs
x = x(:);
x = x(~isnan(x));

% If X is empty, return NaNs
if isempty(x)
    y = NaN(1,class(x));
else
    x = sort(x);
    n = numel(x);    
    if isequal(p,50) % make the median fast
        if rem(n,2) % n is odd
            y = x((n+1)/2);
        else        % n is even
            y = (x(n/2) + x(n/2+1))/2;
        end
    else
        % Compute interpolation
        r = (p/100)*n;
        k = floor(r+0.5);
        if k < 1
            y = x(1);
        elseif k>=n
            y = x(n);
        else
            r = r-k;
            y = (0.5-r).*x(k)+(0.5+r).*x(k+1);
        end
    end

end
