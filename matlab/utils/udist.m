function d2 = udist(u1,u2,lenscale,optimState,dontsum)
%UDIST Vector distance in grid space (needs to be improved, see sq_dist)

if nargin < 5 || isempty(dontsum); dontsum = 0; end

f = optimState.periodicvars;

if any(f)
    D = bsxfun(@minus, u1, u2);
    per = (optimState.UB - optimState.LB)./optimState.scale;
    for d = find(f)
        D(:,d,:) = min(abs(D(:,d,:)), per(d) - abs(D(:,d,:)));
    end    
else
    D = bsxfun(@minus, u1, u2);
end

d2 = bsxfun(@rdivide, D, lenscale).^2;
if ~dontsum; d2 = sum(d2, 2); end

end