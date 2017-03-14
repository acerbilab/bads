function xs = searchGauss(x,gpstruct,LB,UB,optimState,options)
%SEARCHGAUSS Multivariate normal random search.

if nargin == 0
    xs = 'gauss';
    return;
end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

nvars = length(x);

% Diagonal covariance matrix for random search
sigma = 0.5*MeshSize*SearchFactor*eye(nvars);

% Random draw from multivariate normal
xs = bsxfun(@plus, x, randn(options.Nsearch,nvars)*sigma);

end