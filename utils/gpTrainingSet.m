function [gpstruct,exitflag] = gpTrainingSet(gpstruct,method,uc,ui,optimState,options,refit_flag)
%GPTRAININGSET Update training input set of Gaussian Process.
%   GPSTRUCT = GPTRAININGSET(GPSTRUCT,'add',UNEW,FVAL,OPTIMSTATE,OPTIONS)
%   adds point UNEW with function value FVAL to the training set of GP 
%   structure GPSTRUCT. If OPTIONS.NoiseObj is on, FVAL(1) is the estimated
%   function value and FVAL(2) the standard deviation of the estimate.
%   OPTIMSTATE is the optimization structure, and OPTIONS the options
%   structure.
%
%   GPSTRUCT = GPTRAININGSET(GPSTRUCT,METHOD,UC,UI,OPTIMSTATE,OPTIONS)
%   updates training set of GP structure GPSTRUCT according to method 
%   METHOD. UC is the current center/incumbent and UI (optional) pivot
%   vectors used to define the training set. OPTIMSTATE is the optimization 
%   structure, and OPTIONS the options structure.
%
%   GPSTRUCT = GPTRAININGSET(GPSTRUCT,METHOD,UC,UI,OPTIMSTATE,OPTIONS,1)
%   also refits the GP hyperparameters on the training set.
%
%   [GPSTRUCT,EXITFLAG] = GPTRAININGSET(...) returns an EXITFLAG that 
%   describes the exit condition of the GP training. A negative EXITFLAG
%   indicates that the training failed (likely due to error in inverting 
%   the GP covariance matrix at some point during optimization).

%   Luigi Acerbi 2017

exitflag = 0;

if nargin < 7 || isempty(refit_flag); refit_flag = false; end

index = 1:optimState.Xmax;
U = optimState.U(index,:);
Y = optimState.Y(index);
if isfield(optimState,'S'); S = optimState.S(index); end
D = numel(uc);
rotategp_flag = isfield(gpstruct,'C') && ~isempty(gpstruct.C);

switch (lower(method))
    
    case 'add'  % Add a single point to the training set
        if options.FitnessShaping
            ystar = gpstruct.nonlinf(ui(1), gpstruct.nonlinmu, gpstruct.deltay);
        else
            ystar = ui(1);
        end
        try
            if options.NoiseObj; ysd = ui(2); else ysd = []; end
            gpstruct.post = update_posterior(gpstruct.hyp(1), gpstruct.mean, ...
                  gpstruct.cov, gpstruct.x, gpstruct.post, uc, ystar);          
            if rotategp_flag; uc = uc*gpstruct.C'; end     
            gpstruct.x = [gpstruct.x; uc];
            
            % Check that the added value is well-defined
            if isfinite(ystar)
                gpstruct.y = [gpstruct.y; ystar];
                if options.NoiseObj; gpstruct.sd = [gpstruct.sd; ysd]; end
                gpstruct.erry = [gpstruct.erry; false];
            else
                [ypenalty,idx] = max(gpstruct.y);
                gpstruct.y = [gpstruct.y; ypenalty];
                if options.NoiseObj; gpstruct.sd = [gpstruct.sd; gpstruct.sd(idx)]; end
                gpstruct.erry = [gpstruct.erry; true];
            end
            
        catch
            % Posterior update failed, point was not added
        end
        
        return;
          
    case 'nearest' % Add nearest neighbors from incumbent        
        ug = uc;
        
        % Distance between vector and set of reference vectors
        dist = zeros(size(U,1),size(uc,1));
        for i = 1:size(ug,1)
            dist(:,i) = udist(U,ug(i,:),gpstruct.lenscale,optimState);
        end
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');

        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*gpstruct.effectiveradius;
        ntrain = min(options.Ndata, sum(distord <= radius^2));
        
        % Minimum number of points to keep
        ntrain = max([options.MinNdata,options.Ndata-options.BufferNdata,ntrain]);
        
        % Up to the maximum number of available points
        ntrain = min(ntrain, optimState.Xmax);
        
        % Take points closest to reference points
        index = 1:ntrain;
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end
        
        
    case 'grid' % Add points closest to coordinate-wise neighbors
                
        % Add coordinate-wise neighbors to the reference set
        vv = eye(D);
        %vv = grammSchmidt(randn(D,1))';
        vv = [vv; -vv];        
        
        if options.HessianUpdate    % Unsupported
            Bnew = vv*optimState.C';
            % Global vector normalization
            M = sqrt(sum(Bnew(1:D,:).*Bnew(1:D,:),2));
            N = exp(log(M) - mean(log(M)))./M;
            Bnew = bsxfun(@times, Bnew, repmat(N,[size(Bnew,1)/D,1]));
            vv = bsxfun(@times,Bnew*optimState.meshsize,optimState.scale);
            % vv = [vv; bsxfun(@times,[eye(D);-eye(D)]*optimState.meshsize,optimState.scale.*gpstruct.pollscale)];
        else
            vv = bsxfun(@times,vv*optimState.meshsize,optimState.scale.*gpstruct.pollscale);
        end

        ug = periodCheck(bsxfun(@plus,uc,vv),optimState.LB,optimState.UB,optimState);
        ug = uCheck(ug,optimState.TolMesh,optimState,1);
        ug = [uc; ug];
                
        % Distance between vector and set of reference vectors
        dist = zeros(size(U,1),size(ug,1));
        for i = 1:size(ug,1)
            dist(:,i) = udist(U,ug(i,:),gpstruct.lenscale,optimState);
        end
        [~,closest] = min(dist,[],1);        
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');

        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*gpstruct.effectiveradius;
        ntrain = min(options.Ndata, sum(distord <= radius^2));
        
        %------------------------------------------------------------------
        % Consider rescaling the radius by optimState.meshsize        
        %------------------------------------------------------------------

        % Minimum number of points to keep
        ntrain = max([options.MinNdata,options.Ndata-options.BufferNdata,ntrain]);
        
        % Up to the maximum number of available points
        ntrain = min(ntrain, optimState.Xmax);
                
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Take points closest to reference points
        index = 1:ntrain;
        
        % Add closest point
        extraidx = find(any(bsxfun(@eq, ord, closest),2));
        index = [index, extraidx'];
        
        % Add safeguarded points
        for d = 1:D
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
            
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end
        
    %----------------------------------------------------------------------    
    case 'covgrid'
        
        if ~options.HessianUpdate
            error('HESSIANUPDATE must be on.');
        end
        
        % GP-based vector scaling
        vv = eye(D);
        vv = [vv; -vv];
        
        Bnew = vv*optimState.C';
        % Global vector normalization
        M = sqrt(sum(Bnew(1:D,:).*Bnew(1:D,:),2));
        N = exp(log(M) - mean(log(M)))./M;
        Bnew = bsxfun(@times, Bnew, repmat(N,[size(Bnew,1)/D,1]));
        vv = bsxfun(@times,Bnew*optimState.meshsize,optimState.scale);

        %xc = origunits(uc,optimState);
        %xg = periodCheck(bsxfun(@plus,xc,vv),optimState.LB,optimState.UB,optimState,0);
        %xg = xCheck(xg,optimState.LB,optimState.UB,optimState.TolMesh,optimState,0);
        %ug = [uc; gridunits(xg,optimState)];

        ug = periodCheck(bsxfun(@plus,uc,vv),optimState.LB,optimState.UB,optimState);
        ug = uCheck(ug,optimState.TolMesh,optimState,1);
        ug = [uc; ug];        
        
        % Distance between vector and set of poll vectors
        dist = zeros(size(U,1),size(ug,1));
        for i = 1:size(ug,1)
            dist(:,i) = ugdist(U,ug(i,:),optimState.B,optimState);
        end
        [~,closest] = min(dist,[],1);        
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');
        
        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius;
        ntrain = max(options.MinNdata,max(options.Ndata-options.BufferNdata, min(options.Ndata, sum(distord <= radius^2))));
        ntrain = min(ntrain, optimState.Xmax);
        % ntrain = max(floor(ntrain/2),min(ntrain, sum(distord <= radius^2)));
        
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Cluster observations
        index = 1:ntrain;
        
        % Add closest point
        extraidx = find(any(bsxfun(@eq, ord, closest),2));
        index = [index, extraidx'];
        
        % Add safeguarded points
        for d = 1:D
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
            
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);           
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end

    case 'nogrid'
                
        ug = uc;
                
        % Distance between vector and set of poll vectors
        dist = zeros(size(U,1),size(ug,1));
        for i = 1:size(ug,1)
            dist(:,i) = udist(U,ug(i,:),gpstruct.lenscale,optimState);
        end
        [~,closest] = min(dist,[],1);        
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');

        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*gpstruct.effectiveradius;
        ntrain = min(options.Ndata, sum(distord <= radius^2));
        
        % Minimum number of points to keep
        ntrain = max([options.MinNdata,options.Ndata-options.BufferNdata,ntrain]);
        
        % Up to the maximum number of available points
        ntrain = min(ntrain, optimState.Xmax);
                
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Cluster observations
        index = 1:ntrain;
                
        % Add safeguarded points
        for d = 1:D
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
            
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end
        
        
    case 'neighborhood'
        
        if ~options.HessianUpdate
            error('HESSIANUPDATE must be on.');
        end
                        
        dist = ugdist(U,uc,optimState.B,optimState);
        [distord,ord] = sort(dist,'ascend');
        
        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*optimState.meshsize;
        ntrain = max(options.MinNdata,min(options.Ndata, sum(distord <= radius^2)));
        ntrain = min(ntrain, optimState.Xmax);
        
        %distord
        %radius
        
        index = 1:ntrain;                    
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);           
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end        
        
    case 'global'

        globalNdata = min(options.globalNdata, optimState.Xmax);
        
        if optimState.Xmax > globalNdata                
            gpstruct.x = [];
            gpstruct.y = [];
            % Periodic variables are not fully supported for clustering
            idx = kmeans(bsxfun(@rdivide,U,gpstruct.lenscale),globalNdata);
            nClusters = max(idx);
            for iCluster = 1:nClusters
                subset = (idx == iCluster);
                xx = U(subset,:);
                yy = Y(subset);
                if 0
                    gpstruct.y(end+1) = mean(yy);
                    gpstruct.x(end+1,:) = mean(xx,1); 
                else
                    if 1
                        [~,pos] = min(yy);
                    else
                        pos = randi(numel(yy));
                    end
                    gpstruct.y(end+1) = yy(pos);
                    gpstruct.x(end+1,:) = xx(pos,:); 
                end
            end
        else
            gpstruct.x = U;
            gpstruct.y = Y;            
        end
        
        gpstruct.y = gpstruct.y(:);
    %----------------------------------------------------------------------
end

% Transformation of objective function
if options.FitnessShaping
    [gpstruct.y,gpstruct.nonlinf,gpstruct.nonlinmu,gpstruct.deltay] = ...
        fitnessTransform(gpstruct.y);
end

% Substitute ill-defined function values with highest value in set
err_index = ~isfinite(gpstruct.y);
if any(err_index)
    [ypenalty,idx1] = max(gpstruct.y(~err_index));
    idx_values = find(~error_index);
    gpstruct.y(err_index) = ypenalty;
    if isfield(optimState,'S'); gpstruct.sd(err_index) = gpstruct.sd(idx_values(idx1)); end
end

% List of 'erroneous' points
gpstruct.erry = err_index;

% Store test points
gpstruct.xi = ui;

% Rotate dataset (unsupported)
if rotategp_flag && ~strcmpi(method,'add') && ~strcmpi(method,'neighborhood')
    if refit_flag
        if isfield(gpstruct,'Cinv')
            Cinvold = gpstruct.Cinv;
        else
            Cinvold = eye(D);
        end
        C = inv(optimState.C)*optimState.meshsize;
        gpstruct.C = C;
        gpstruct.Cinv = optimState.C/optimState.meshsize;
        idx = gpstruct.ncovoffset+(1:D);
        covlen = exp(gpstruct.hyp.cov(idx));
        covlen = log(gpstruct.C*Cinvold*covlen) % This needs to be fixed
        for i = 1:D
            covlen(i) = min(max(covlen(i), gpstruct.bounds.cov{i}(1)),gpstruct.bounds.cov{i}(2));
        end
        gpstruct.hyp.cov(idx) = covlen;
        
        % C
    end
    gpstruct.x = gpstruct.x*gpstruct.C';
end

% Update GP definitions
gplik = [];
gpstruct.x0 = [];
gpstruct = feval(options.gpdefFcn{:},D,gplik,optimState,options,gpstruct);

% Re-fit Gaussian process (optimize or sample -- only optimization supported)
if refit_flag
    
    [gpstruct,exitflag] = gpfit(gpstruct,options.gpSamples,options);
    Nsamples = numel(gpstruct.hyp);
    
    % Gaussian process length scale
    if gpstruct.ncovlen > 1
        gpstruct.lenscale = zeros(1,D);
        for i = 1:Nsamples
            gpstruct.lenscale = gpstruct.lenscale + gpstruct.hypweight(i)*exp(gpstruct.hyp(i).cov(gpstruct.ncovoffset+(1:D)))';
        end
    else
        gpstruct.lenscale = 1;
    end
    
    % GP-based geometric length scale
    ll = zeros(D,Nsamples);
    for i = 1:Nsamples
        ll(:,i) = options.gpRescalePoll*gpstruct.hyp(i).cov((1:D)+gpstruct.ncovoffset);
    end    
    ll = exp(sum(bsxfun(@times, gpstruct.hypweight, ll - mean(ll(:))),2))';
    
    % Take bounded limits
    ub_bounded = optimState.UB;
    ub_bounded(~isfinite(ub_bounded)) = optimState.PUB(~isfinite(ub_bounded));
    lb_bounded = optimState.LB;
    lb_bounded(~isfinite(lb_bounded)) = optimState.PLB(~isfinite(lb_bounded));    
    
    ll = min(max(ll, optimState.searchmeshsize), (ub_bounded-lb_bounded)./optimState.scale); % Perhaps this should just be PUB - PLB?
    gpstruct.pollscale = ll;
    
    % GP effective covariance length scale radius
    if options.UseEffectiveRadius
        switch lower(gpstruct.covtype)
            case 'rq'
                alpha = zeros(1,Nsamples);
                for i = 1:Nsamples; alpha(i) = exp(gpstruct.hyp(i).cov(end)); end
                alpha = sum(gpstruct.hypweight.*alpha);
                gpstruct.effectiveradius = sqrt(alpha*(exp(1/alpha)-1));                
            case 'matern1'
                gpstruct.effectiveradius = 1/sqrt(2);
            case 'matern3'
                % gpstruct.effectiveradius = fzero(@(x)(1+sqrt(3)*x)*exp(-sqrt(3)*x)-exp(-1),1)/sqrt(2);
                gpstruct.effectiveradius = 0.876179713323485;
            case 'matern5'
                % gpstruct.effectiveradius = fzero(@(x)(1+sqrt(5)*x+5/3*x^2)*exp(-sqrt(5)*x)-exp(-1),1)/sqrt(2);
                gpstruct.effectiveradius = 0.918524648109253;
            otherwise
                gpstruct.effectiveradius = 1;
        end
    end
    
    % Gaussian process signal variability
    gpstruct.sf = 0;
    for i = 1:Nsamples
        gpstruct.sf = gpstruct.sf + gpstruct.hypweight(i)*exp(gpstruct.hyp(i).cov(gpstruct.ncovoffset+D+1));
    end
    
    % [std(gpstruct.y) gpstruct.sf gpstruct.hyp.lik(1) exp(gpstruct.hyp.lik(2:end)')]
end

try
    % Recompute posterior
    Nsamples = numel(gpstruct.hyp);    
    [~,~,~,~,~,gpstruct.post] = mygp(gpstruct.hyp(1),gpstruct.inf,gpstruct.mean,gpstruct.cov,gpstruct.lik,gpstruct.x,gpstruct.y,uc(1,:));    
    for i = 2:Nsamples
        [~,~,~,~,~,gpstruct.post(i)] = mygp(gpstruct.hyp(i),gpstruct.inf,gpstruct.mean,gpstruct.cov,gpstruct.lik,gpstruct.x,gpstruct.y,uc(1,:));
    end
catch
    % Posterior update failed
    gpstruct.post = [];
    exitflag = -2;
end

% gpstruct.hyp.cov(:)'
    
end

%--------------------------------------------------------------------------
function [gpstruct,exitflag] = gpfit(gpstruct,Nsamples,options)
%GPFIT Fit Gaussian Process hyper-parameters (optimize or sample).

if isfield(gpstruct,'bounds') && ~isempty(gpstruct.bounds)
    bounds = unwrap2vec(gpstruct.bounds);
    lb = bounds(1:2:end-1);
    ub = bounds(2:2:end);
else
    lb = -Inf;
    ub = Inf;
end

% First find MAP

% Initial point #1 is old hyperparameter value (randomly picked)
hyp0(1) = gpstruct.hyp(randi(numel(gpstruct.hyp)));

% Initial point #2 is avg of random draw from prior and #1

% Check for possible high-noise mode
if numel(options.NoiseSize) == 1 || ~isfinite(options.NoiseSize(2)); noise = 1; 
else noise = options.NoiseSize(2); end    
highNoise = hyp0(1).lik(1) > (log(options.NoiseSize(1)) + 2*noise);

% Check for mean stuck below minimum
lowMean = hyp0(1).mean(1) < min(gpstruct.y);

% Conditions for performing a second fit
secondfit = options.DoubleRefit || highNoise || lowMean;

if secondfit
    hrnd = gppriorrnd(gpstruct.prior,gpstruct.hyp(1));
    hrnd = 0.5*(unwrap2vec(hrnd) + unwrap2vec(gpstruct.hyp(1)));
    if highNoise; hrnd(end-1) = randn()-2; end   % Retry with low noise magnitude
    if lowMean; hrnd(end) = median(gpstruct.y); end % Retry with mean from median
    hyp0(2) = rewrap(gpstruct.hyp(1),min(max(hrnd,lb),ub));
end
optoptions = optimset('TolFun',0.1,'TolX',1e-4,'MaxFunEval',150);
% optoptions.Hessian = 'user-supplied';

% Remove undefined points
if isfield(gpstruct,'erry') && sum(gpstruct.erry) > 0 && 0
    gpopt = gpstruct;
    gpopt.x(gpstruct.erry,:) = [];
    gpopt.y(gpstruct.erry) = [];
    if isfield(gpopt.sd); gpopt.sd(gpstruct.erry) = []; end
else
    gpopt = gpstruct;
end

[hyp,exitflag] = gpHyperOptimize(hyp0,gpopt,options.OptimToolbox,optoptions,options.NoiseNudge,options.RemovePointsAfterTries);     
hypw = 1;

% If using multiple samples, do SVGD from neighborhood of the MAP
if Nsamples > 1
    [hyp,hypw] = gpHyperSVGD(hyp,gpstruct,options);    
end

gpstruct.hyp = hyp;
gpstruct.hypweight = hypw;         % Update samples weigths

% Laplace approximation at MAP solution
% if Nsamples == 0 && gpstruct.marginalize
%     try
%         [~,~,~,~,~,gpstruct.hypHessian] = feval(gpstruct.inf{:},gpstruct.hyp,gpstruct.mean,gpstruct.cov,gpstruct.lik,gpstruct.x,gpstruct.y);
%     catch
%         gpstruct.hypHessian = [];
%     end
% end

%[exp(gpstruct.hyp.cov(:))', exp(gpstruct.hyp.lik)]
%gpstruct.hyp.mean

end