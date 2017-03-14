function xs = searchWCM(sumrule,x,gpstruct,LB,UB,optimState,options)
%SEARCHWCM Weighted covariance matrix search step (inspired by CMA-ES).

if nargin < 2
    xs = 'wcm';
    return;
end

% SUMRULE not specified, all arguments shifted by one
if nargin < 7
    options = optimState;
    optimState = UB;
    UB = LB;
    LB = gpstruct;
    gpstruct = x;
    x = sumrule;
    sumrule = [];    
end

% By default renormalize by expected vector magnitude
if isempty(sumrule); sumrule = 1; end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

X = gpstruct.x;
Y = gpstruct.y;

nvars = length(x);

% Small jitter added to each direction
% jit = (options.PollMeshMultiplier^-options.SearchGridNumber)*MeshSize;
jit = MeshSize;

% Compute vector weights
mu = 0.5*size(X,1);
weights = zeros(1,1,floor(mu));
weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
weights = weights./sum(weights);

% Compute top vectors
[~,index] = sort(Y,'ascend');

% Compute weighted covariance matrix wrt X0
topx = bsxfun(@minus,X(index(1:floor(mu)),:),x);
C = sum(bsxfun(@times,weights,topx'*topx),3);

% Rescale covariance matrix according to mean vector length
[E,lambda] = eig(C);
% [sqrt(diag(lambda))',jit]
lambda = diag(lambda) + jit^2;
if sumrule; lambda = lambda/sum(lambda); else lambda = lambda/max(lambda); end

% Square root of covariance matrix
sigma = diag(sqrt(lambda))*E';

% Rescale by current scale (reduced)
sigma = MeshSize*SearchFactor*sigma;

% Random draw from multivariate normal
xs = bsxfun(@plus, x, randn(options.Nsearch,nvars)*sigma);

% xs = [xs; bsxfun(@plus, x, randn(options.Nsearch,nvars)*eye(nvars)/sqrt(nvars)*0.5*MeshSize*SearchFactor)];

if 0    % Search only a random subset of coordinates?
    ndim = ceil(sqrt(nvars));
    index = randperm(nvars);
    R = zeros(options.Nsearch,nvars);
    R(:,index(1:ndim)) = randn(options.Nsearch,ndim);
    xs = bsxfun(@plus, x, R*sigma);
end

end