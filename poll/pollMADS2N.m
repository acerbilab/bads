function Bnew = pollMADS2N(B,x,gpstruct,LB,UB,optimState,options)
%POLLMADS2N Poll 2N random basis vectors (mesh-adaptive direct search).

nvars = numel(x);

if isempty(B)
    nmax = max(1, round(optimState.searchmeshsize/optimState.meshsize));
    
    if nmax > 0
        D = tril(randi(nmax*2-1,nvars) - nmax,-1);
    else
        D = zeros(nvars);
    end
    D = D + diag(nmax*2*(randi(2,[1, nvars])-1.5));
    
    % Random permutation of rows and columns and then transpose
    D = D(randperm(nvars),randperm(nvars))';
    
    % D(abs(D) < nmax/2) = 0;
        
    % sqrt(sum(D.^2,2))*optimState.meshsize
    
    % Counteract subsequent multiplication by pollscale
    D = bsxfun(@rdivide, D, gpstruct.pollscale);
        
    % Axis-oriented poll vector set
    Bnew = [D; -D];
    
    
else
    Bnew = [];
end

end