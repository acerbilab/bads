function [us,optimState] = searchHedge(searchFcns,u,gpstruct,LB,UB,optimState,options)
%SEARCHWCM Weighted covariance matrix search step (inspired by CMA-ES).

if nargin < 2 || isempty(u) || isempty(gpstruct)
    try
        us = feval(searchFcns{optimState.hedge.chosen}{:});
    catch
        us = 'hedge';
    end
    return;
end

% SUMRULE not specified, all arguments shifted by one
if nargin < 7
    options = optimState;
    optimState = UB;
    UB = LB;
    LB = gpstruct;
    gpstruct = u;
    u = sumrule;
    sumrule = [];    
end


hedge = optimState.hedge;

nh = numel(searchFcns);

% Initialize hedge struct
if isempty(hedge)
    hedge.g = zeros(1,nh);
    hedge.g(1) = 10;
    for i = 1:nh; hedge.str{i} = feval(searchFcns{i}{:}); end    
    hedge.n = numel(hedge.g);
    hedge.count = 0;
    hedge.gamma = options.HedgeGamma;
    hedge.beta = options.HedgeBeta;
    hedge.decay = options.HedgeDecay;
    % hedge.string = {'mpi','mei'};
end

hedge.count = hedge.count + 1;
% gammaHedge = min(1/nHedge, sqrt(log(nHedge)/(nHedge*hedge.count)));

hedge.p = exp(hedge.beta*(hedge.g - max(hedge.g)))./sum(exp(hedge.beta*(hedge.g - max(hedge.g))));
hedge.p = hedge.p*(1-hedge.n*hedge.gamma) + hedge.gamma;
% hedge.p
hedge.chosen = find(rand() < cumsum(hedge.p),1);
if isempty(hedge.chosen)
    hedge.chosen = randi(nh);
    if options.gpWarnings
        warning('bads:searchHedgeFail', ['Cannot determine best search function in SEARCHHEDGE (P=' numarray2str(hedge.p) '). Attempting to continue.']);
    end
end
if hedge.gamma == 0
    hedge.phat = ones(size(hedge.g));
else
    hedge.phat = Inf(size(hedge.p));
    hedge.phat(hedge.chosen) = hedge.p(hedge.chosen);
end

us = feval(searchFcns{hedge.chosen}{:}, ...
    u, ...
    gpstruct, ...
    LB, ...
    UB, ...
    optimState, ...
    options);

optimState.hedge = hedge;



end