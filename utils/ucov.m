function C = ucov(U,u0,w,optimState)
%UCOV Weighted covariance matrix of vectors U with respect to U0
        
if isempty(w)
    weights = 1;
else
    weights(1,1,:) = w(:);
end

% Shift periodic variables
f = optimState.periodicvars;
per = (optimState.UB - optimState.LB)./optimState.scale;
for d = find(f)
    U(:,d) = mod(U(:,d) - u0(d) + 0.5*per(d),per(d)) - 0.5*per(d);
    u0(d) = 0;
end

ushift = bsxfun(@minus,U,u0);
C = sum(bsxfun(@times,weights,ushift'*ushift),3);

end