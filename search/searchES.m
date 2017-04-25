function [us,optimState] = searchES(method,sumrule,u,gpstruct,LB,UB,optimState,options)
%SEARCHES Evolution strategies search.

if nargin < 3 || isempty(u) || isempty(gpstruct)
    switch method
        case 1; us = 'ES-wcm';
        case 2; us = 'ES-ell';
        case 3; us = 'ES-eye';
        case 4; us = 'ES-cov';
        case 5; us = 'ES-cma+';
        otherwise; us = 'ES';
    end
    return;
end

% SUMRULE not specified, all arguments shifted by one
if nargin < 8
    options = optimState;
    optimState = UB;
    UB = LB;
    LB = gpstruct;
    gpstruct = u;
    u = sumrule;
    sumrule = [];    
end

% By default renormalize by expected vector magnitude
if isempty(sumrule); sumrule = 1; end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;
rotategp_flag = isfield(gpstruct,'C') && ~isempty(gpstruct.C);

U = gpstruct.x;
Y = gpstruct.y;

nvars = length(u);

switch method

    case {1,5}
        % Small jitter added to each direction
        if method == 1
            jit = MeshSize;
            active_flag = false;
            frac = 0.5;
        else
            % jit = sqrt(optimState.searchmeshsize*MeshSize);
            jit = optimState.searchmeshsize;
            active_flag = true;
            frac = 0.25;
        end

        % Compute vector weights
        mu = frac*size(U,1);
        weights = zeros(1,1,floor(mu));
        weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
        weights = weights./sum(weights);

        % Compute best vectors
        [~,index] = sort(Y,'ascend');
        Ubest = U(index(1:floor(mu)),:);        
        if active_flag; Uworst = U(index(end:-1:end-floor(mu)+1),:); end

        % Compute weighted covariance matrix wrt u0
        C = ucov(Ubest,u,weights,optimState);
        if active_flag
            negC = ucov(Uworst,u,weights,optimState);
            negmueff = sum(1./weights.^2);            
            negcov = 0.25 * negmueff / ((nvars+2)^1.5 + 2*negmueff);
            % [E,lambda] = eig(negC);
            % D = lambda;
            % C = C*((C + inv(negC))\inv(negC));
            C = C - negcov*negC;
        end

        % Rescale covariance matrix according to mean vector length
        [E,lambda] = eig(C);
        % [sqrt(diag(lambda))',jit]
        lambda = diag(max(0,lambda)) + jit^2;
        if sumrule; lambda = lambda/sum(lambda); else lambda = lambda/max(lambda); end

        % Square root of covariance matrix
        sqrtsigma = diag(sqrt(lambda))*E';

    case 2
        % Recaled length scale based on gp
        rescaledLenscale = gpstruct.pollscale;
        rescaledLenscale = rescaledLenscale/sqrt(sum(rescaledLenscale.^2));
        if rotategp_flag
            sqrtsigma = gpstruct.Cinv'*diag(rescaledLenscale);
        else
            sqrtsigma = diag(rescaledLenscale);
        end

    case 3  
        sqrtsigma = eye(nvars)/sqrt(nvars);
        
    case 4
        sqrtsigma = optimState.C';        
end


% Rescale by current scale
sqrtsigma = MeshSize*SearchFactor*sqrtsigma;    

N = optimState.es.mu;

vec = {[-1,0],[-1,0],[-1,0],[-1,0],[-1,0]};
w = options.PollMeshMultiplier.^vec{method};
ns = diff(round(linspace(0,N,numel(w)+1)));

v = [];
for i = 1:numel(w); v = [v; w(i)*ones(ns(i),1)]; end

% Initial population
unew = bsxfun(@plus, u, bsxfun(@times, v, randn(N,nvars)*sqrtsigma));

scale = options.ESstart;

us = []; zold = [];
lambda = optimState.es.lambda;

% Loop over evolutionary strategies iterations
for i = 1:options.Nsearchiter
    
    % Enforce periodicity
    unew = periodCheck(unew,LB,UB,optimState);
    
    % Force candidate points on search grid
    unew = force2grid(unew,optimState);

    % Remove already evaluated or unfeasible points from search set
    unew = uCheck(unew,optimState.TolMesh,optimState,1);
    
    z = [];
    try
        % tic
        if options.AcqHedge
            error('Hedge not supported here.');
            %[optimState.hedge,acqIndex,ymu(:,idx),ys(:,idx)] = ...
            %    acqPortfolio('acq',optimState.hedge,usearch(idx,:),ftarget,fstarget,gpstruct,optimState,options,SufficientImprovement);
            %index = acqIndex(optimState.hedge.chosen);
            %z = 1;
        else
            % Batch evaluation of acquisition function on search set
            [z,~,~,~,fmu,fs] = feval(options.SearchAcqFcn{:},unew,optimState.ftarget,gpstruct,optimState,0);
        end
    catch
        % Something went wrong...
        fmu = optimState.ftargetmu;        
    end
    
    AcqFun = options.SearchAcqFcn{1};
    if isa(AcqFun, 'function_handle'); AcqFun = func2str(AcqFun); end    
    if ischar(AcqFun) && any(strcmp(AcqFun, {'acqNegEIMin','acqNegPIMin'}))
        eta = options.SearchAcqFcn{2};
        if isempty(eta); eta = 0.1; end
        [~, idx] = min([optimState.ftargetmu - eta*optimState.ftargets; fmu(:) - eta*fs(:)]);
        if idx > 1 && ~isempty(zold)  % New minimum found, recompute zold
            optimState.ftargetmu = fmu(idx-1);
            optimState.ftargets = fs(idx-1);
            zold = feval(options.SearchAcqFcn{:},us,optimState.ftarget,gpstruct,optimState,0);
        end
    end
    
    % Something went wrong, random search
    if isempty(z); z = rand(size(unew,1),1); end
    
    nold = size(us,1);
    us = [us; unew];
    z = [zold(:); z(:)];

    N = min(size(us,1),optimState.es.lambda);
    
    % Order candidates and select
    [z,index] = sort(z);
    % [size(xnew,1),nold]
    ntest = min(size(unew,1),nold);
    nnew = sum(index(1:ntest) > nold);
    zold = z(1:N);
    us = us(index(1:N),:);
        
    if ~isempty(zold)
        zlist(i) = zold(1); % Store best Z value at i-th iteration
    else
        zlist(i) = NaN;
    end
    
    if i < options.Nsearchiter
        % Update scale parameter
        frac = nnew/ntest;
        if i > 1; scale = scale*exp(options.ESbeta*(frac-0.2)); end
        % [i frac scale]
        
        % Reproduce
        
        es = ESupdate(size(us,1),lambda,optimState.es.iter);        
        ll = min(lambda, size(us,1));
        % xnew = xs(optimState.es.selectmask(1:ll),:) + randn(ll, nvars)*sigma*scale;
        unew = us(es.selectmask(1:ll),:) + randn(ll, nvars)*sqrtsigma*scale;
    end
    
end

% zlist

% Return best
if ~isempty(us); us = us(1,:); end

end