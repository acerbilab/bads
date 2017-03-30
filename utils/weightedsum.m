function z = weightedsum(w,y,dim)
%WEIGHTEDSUM Perform weighted sum.

if nargin < 3 || isempty(dim); dim = 1; end

switch dim
    case 1
        w = w(:);
        z = sum(bsxfun(@times,w(~isnan(w)),y(~isnan(w),:)),1);
    case 2
        w = w(:)';
        z = sum(bsxfun(@times,w(~isnan(w)),y(:,~isnan(w))),2);
        
    otherwise
        error('Summation dimension DIM needs to be 1 or 2.');
end