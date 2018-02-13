function hypr = gppriorrnd(prior,hyp)
%GPPRIORRND Random draw from GP prior.
%
%  HYPR = GPPRIORRND(PRIOR,HYP) draws a random set of hyperparameters 
%  from GP prior PRIOR, returned as a hyperparameter struct arranged
%  according to struct HYP. 
%
%  A hyperparameter with an empty prior, @priorDelta (or 'priorDelta'),
%  or @priorClamped (or 'priorClamped') is taken unchanged from HYP.
%
%  See also INFPRIOR, PRIORDISTRIBUTION, USAGEPRIOR.
%
%  Luigi Acerbi, 2015-Oct-03 
%  Modified from infPrior.m by Hannes Nickisch and Roman Garnett.
%
% See also INFMETHODS.M, USAGEPRIOR.M, PRIORDISTRIBUTIONS.M.

if isempty(prior); return; end

hypr = hyp;
nam = fieldnames(prior);
num = zeros(numel(unwrap2vec(hyp)),1);       % number of priors per hyperparameter

for i=1:numel(nam)     % iterate over kinds of hyperparameters cov/lik/mean/xu
    ni = nam{i};                                         % name of the ith field
    
    % first catch multivariate priors
    if strcmp(ni,'multi')
        p = prior.(ni);
        for j=1:numel(p)  % iterate over individual muti-hyperparameter components
            pj = p{j};
            if ~isempty(pj)          % only proceed if a nonempty prior is specified
                idx = pj{end}; pj = pj(1:end-1);            % grab index from cell end
                pjstr = pj{1}; if ~ischar(pjstr), pjstr = func2str(pjstr); end
                if numel(pjstr)<5 || ~strcmp(pjstr(end-4:end),'Multi')
                    error('multivariate hyperpriors are called <Name>Multi')
                end
                if isstruct(idx)                           % massage idx into a vector
                    idxs = rewrap(hypr,zeros(size(num)));              % structured index
                    for nj = fieldnames(idx), idxs.(nj{1})( idx.(nj{1}) ) = 1; end
                    idx = unwrap2vec(idxs)>0;                   % linearise structured index
                else
                    idxz = zeros(size(num)); idxz(idx) = 1; idx = idxz>0; % binary index
                end
                if sum(idx)<=1, error('multivariate priors need >1 hyperparam'), end
                num(idx) = num(idx)+1;                             % inc prior counter
                hypu = unwrap2vec(hypr);
                if strncmp(pjstr,'priorClamped',12) || ...
                    strncmp(pjstr,'priorDelta',10)
                        r(idx) = hypu(idx);
                else
                    r(idx) = feval(pj{:});    % evaluate prior distribution
                end
            else
              error('multivariate priors should be non empty')
            end
        end
        continue                                       % jump to univariate priors
    end
    
    if ~isfield(hyp,ni), error(['unknown prior field ',ni,' in hyp']), end
    
    p = prior.(ni);
    if numel(p)~=numel(hyp.(ni)), error(['bad hyp/prior field length ',ni]), end
    
    for j=1:numel(p)         % iterate over individual hyperparameter components
        pj = p{j}; if ~iscell(pj) && ~isempty(pj), pj = {pj}; end   % enforce cell
        if ~isempty(pj)            % only proceed if a nonempty prior is specified
            num = rewrap(hypr,num); num.(ni)(j) = num.(ni)(j)+1;  % inc prior counter
            num = unwrap2vec(num);
            pj1str = pj{1}; if ~ischar(pj1str), pj1str = func2str(pj1str); end
            if strncmp(pj1str,'priorClamped',12) || strncmp(pj1str,'priorDelta',10)
                hypr.(ni)(j) = hyp.(ni)(j);     % copy value from base          
            else
                hypr.(ni)(j) = feval(pj{:});    % random sample from prior distribution
            end
        end
    end
end

% Check for hypers with more than a single prior
if any(num>1)                 
    num = rewrap(hypr,num); nam = fieldnames(num);
    s = '';
    for i=1:numel(nam)
        idx = find(num.(nam{i})>1);
        for j=1:numel(idx)
            s = [s,sprintf('hyp.%s(%d) ',nam{i},idx(j))];
        end
    end
    error(['More than 1 prior specified for ',s(1:end-1),'.'])  
end
