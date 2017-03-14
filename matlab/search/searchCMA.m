function [xs,optimState] = searchCMA(sumrule,x,gpstruct,LB,UB,optimState,options)
%SEARCHCMA Covariance matrix adaptation search step (inspired by CMA-ES).

if nargin < 2 || isempty(x) || isempty(gpstruct)
    xs = 'cma';
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

nvars = length(x);

sigma = optimState.C';
sigma = MeshSize*SearchFactor*sigma;    % Rescale by current scale

N = options.Nsearch;
xs = bsxfun(@plus, x, randn(N,nvars)*sigma);

end