function d2 = ugdist(u1,u2,B,optimState,dontsum)
%UGDIST Vector distance in grid space with generic bilinear form B

if nargin < 5 || isempty(dontsum); dontsum = 0; end

f = optimState.periodicvars;

if any(f)
    error('Periodic variables not supported yet.');
    D = bsxfun(@minus, u1, u2);
    per = (optimState.UB - optimState.LB)./optimState.scale;
    for d = find(f)
        D(:,d,:) = min(abs(D(:,d,:)), per(d) - abs(D(:,d,:)));
    end    
    d2 = bsxfun(@rdivide, D, lenscale).^2;
else
    if size(u2,1) == 1; u2 = repmat(u2,[size(u1,1),1]); end    
    d2 = diag(u1*B*u1' + u2*B*u2' - 2*u1*B*u2');
end

if ~dontsum; d2 = sum(d2, 2); end

end