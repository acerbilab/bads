function X = xCheck(X,LB,UB,TolMesh,optimState,xgridscaled)
%XCHECK Remove duplicates, vectors outside bounds or previously evaluated

% If X is already rescaled, rescale also LB and UB
if xgridscaled
    UB = gridunits(UB,optimState);
    LB = gridunits(LB,optimState);
end

f = logical(sum(bsxfun(@gt,X,UB) | bsxfun(@lt,X,LB),2));

if 0
    % Project vectors outside bounds on bounds
    if 0
        if xgridscaled
            x = gridunits(optimState.x,optimState);
        else
            x = optimState.x;
        end

        if any(f)
            dv = bsxfun(@minus, X(f,:), optimState.x);
            nf = sqrt(sum(dv.*dv,2));  % normalization
            du = bsxfun(@rdivide, dv, nf);

            up = bsxfun(@rdivide, UB - x, du);
            down = bsxfun(@rdivide, x - LB, -du);

            r = up.*(du >= 0) + down.*(du < 0);
            rmin = max(0, min(r,[],2) - optimState.searchmeshsize/2);

            Xf = bsxfun(@plus, x, bsxfun(@times, rmin, du));
            Xf = optimState.searchmeshsize.*round(bsxfun(@rdivide, Xf, optimState.searchmeshsize));

            X(f,:) = Xf;
            f = logical(sum(bsxfun(@gt,X,UB) | bsxfun(@lt,X,LB),2));
        end
    else
        tol = TolMesh;
        X = bsxfun(@max, bsxfun(@min, X, UB - 0.5*tol), LB + 0.5*tol);
        X = tol.*round(bsxfun(@rdivide, X, tol));
    end
end

X(f,:) = [];

% Remove duplicate vectors
X = unique(X,'rows');

% Remove previously evaluted vectors (within TolMesh)
if ~isempty(X)
    TolMesh = TolMesh/2;
    % TolMesh = exp(0.5*(log(TolMesh) + log(optimState.meshsize)))/2;
    % Convert everything to normalized grid coordinates for comparison
    if xgridscaled
        u1 = round(X./TolMesh);        
    else
        u1 = round(gridunits(X,optimState)./TolMesh);
    end
    U = optimState.U(1:optimState.Xmax,:);
    u2 = round(U./TolMesh);
    [~,idx] = setdiff(u1,u2,'rows');
    X = X(idx,:);
end

end