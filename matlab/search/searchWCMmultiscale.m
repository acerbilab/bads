function [xs,optimState] = searchWCMmultiscale(method,sumrule,x,gpstruct,LB,UB,optimState,options)
%SEARCHWCM Weighted covariance matrix search step (inspired by CMA-ES).

if nargin < 3 || isempty(x) || isempty(gpstruct)
    xs = ['wcm+' num2str(method)];
    return;
end

% SUMRULE not specified, all arguments shifted by one
if nargin < 8
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

if ~isfield(optimState,'searchsigma') || isempty(optimState.searchsigma) || 1

    switch method

        case 1
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
            xmu = x;
            topx = bsxfun(@minus,X(index(1:floor(mu)),:),xmu);
            C = sum(bsxfun(@times,weights,topx'*topx),3);

            % Rescale covariance matrix according to mean vector length
            [E,lambda] = eig(C);
            % [sqrt(diag(lambda))',jit]
            lambda = diag(lambda) + jit^2;
            if sumrule; lambda = lambda/sum(lambda); else lambda = lambda/max(lambda); end

            % Square root of covariance matrix
            sigma = diag(sqrt(lambda))*E';

        case 2
            % Recaled length scale based on gp
            rescaledLenscale = gpstruct.pollscale;
            rescaledLenscale = rescaledLenscale/sqrt(sum(rescaledLenscale.^2));
            sigma = diag(rescaledLenscale);

        case 3  
            sigma = eye(nvars)/sqrt(nvars);
            
        case 4
            sigma = optimState.C';
            
            
    end

    optimState.searchsigma = sigma;    
else
    sigma = optimState.searchsigma;
    optimState.searchsigma = [];    
end

% Rescale by current scale
sigma = MeshSize*SearchFactor*sigma;    

N = round(options.Nsearch);

vec = {[-1,0],[-1,0],[-1,0]};
w = options.PollMeshMultiplier.^vec{method};
ns = diff(round(linspace(0,N,numel(w)+1)));

v = [];
for i = 1:numel(w); v = [v; w(i)*ones(ns(i),1)]; end

xs = bsxfun(@plus, x, bsxfun(@times, v, randn(N,nvars)*sigma));

end