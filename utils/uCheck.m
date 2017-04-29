function U = uCheck(U,TolMesh,optimState,proj)
%UCHECK Remove duplicates, vectors outside bounds or previously evaluated

if nargin < 4 || isempty(proj); proj = 0; end

if proj
    % Project vectors outside bounds on search mesh points closest to bounds
    U = bsxfun(@max, bsxfun(@min, U, optimState.UBsearch), optimState.LBsearch);
else
    idx = any(bsxfun(@gt,U,optimState.UB) | bsxfun(@lt,U,optimState.LB),2);
    U(idx,:) = [];
end

% Remove duplicate vectors
U = unique(U,'rows');

% Remove previously evaluted vectors (within TolMesh)
if ~isempty(U)
    TolMesh = TolMesh/2;
    % TolMesh = exp(0.5*(log(TolMesh) + log(optimState.meshsize)))/2;
    % Convert everything to normalized grid coordinates for comparison
    u1 = round(U./TolMesh);        
    U2 = optimState.U(1:optimState.Xmax,:);
    u2 = round(U2./TolMesh);
    [~,idx] = setdiff(u1,u2,'rows');
    U = U(idx,:);
end

% Evaluate non-bound constraints
if ~isempty(optimState.nonbcon)
    X = origunits(U,optimState);    % Convert back to original space
    C = optimState.nonbcon(X);      % Evaluate constraints
    idx = (C <= 0);                 % Keep points that satisfy constraints
    U = U(idx,:);
end

end