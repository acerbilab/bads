function gpstruct = gpdefBads(covtype,covextras,D,gplik,optimState,options,gpstruct)
%GPDEFBADS Define Gaussian Process (GP) structure for BADS.

if nargin < 7; gpstruct = []; end

MeshSize = optimState.meshsize;
TolFun = options.TolFun;
TolMesh = optimState.TolMesh;

gaussPriorFunc = @priorGauss;

covard = covextras(1);   % First covariance flag (ARD covariance)

% Initialize new GP struct
if isempty(gpstruct)

    %% GP covariance function
    gpstruct.covtype = lower(covtype);
        
    switch lower(covtype)
        case 'se'
            if covard;  gpstruct.cov = {@covSEard_fast};
            else        gpstruct.cov = {@covSEiso};         end
        case 'matern1'
            if covard;  gpstruct.cov = {@covMaternard_fast,1};
            else        gpstruct.cov = {@covMaterniso,1};   end
        case 'matern3'
            if covard;  gpstruct.cov = {@covMaternard_fast,3};
            else        gpstruct.cov = {@covMaterniso,3};   end
        case 'matern5'
            if covard;  gpstruct.cov = {@covMaternard_fast,5};
            else        gpstruct.cov = {@covMaterniso,5};   end
        case 'pp0'
            if covard;  gpstruct.cov = {@covPPard,0};
            else        gpstruct.cov = {@covPPiso,0};   end            
        case 'pp1'
            if covard;  gpstruct.cov = {@covPPard,1};
            else        gpstruct.cov = {@covPPiso,1};   end            
        case 'pp2'
            if covard;  gpstruct.cov = {@covPPard,2};
            else        gpstruct.cov = {@covPPiso,2};   end            
        case 'pp3'
            if covard;  gpstruct.cov = {@covPPard,3};
            else        gpstruct.cov = {@covPPiso,3};   end            
        case 'rq'
            if covard;  gpstruct.cov = {@covRQard_fast};
            else        gpstruct.cov = {@covRQiso};   end
        case 'lin'
            if covard;  gpstruct.cov = {@covSum,{@covConst,@covLINard}};
            else        gpstruct.cov = {@covRQiso};   end            
        otherwise
            error(['Unknown covariance type ''' covtype ''' in GP definition.']);
    end
        
    % Total # of covariance hyperparameters
    ncov = eval(feval(gpstruct.cov{:}));
    gpstruct.hyp.cov = zeros(ncov,1);
    
    % # of length scale hyperparameters
    if strcmpi(covtype,'rq'); ncovlen = ncov - 2; 
    else ncovlen = ncov - 1; end
    gpstruct.ncovlen = ncovlen;
    
    gpstruct.prior.cov = [];
    gpstruct.bounds.cov = [];
    
    % Prior and bounds for covariance periodic length
    if any(optimState.periodicvars)
        per = log((optimState.UB - optimState.LB)./optimState.scale);
        
        % Infinite period for non-periodic dimensions
        per(~optimState.periodicvars) = Inf;
                
        % Add period to covariance hyperparameters
        gpstruct.hyp.cov = [per(:); gpstruct.hyp.cov];
        
        for i = 1:gpstruct.ncovlen
            gpstruct.prior.cov{end+1} = {@priorDelta};  % Period is fixed
            gpstruct.bounds.cov{end+1} = [-Inf; Inf];
        end
        gpstruct.ncovoffset = gpstruct.ncovlen;
        
        if ncovlen == 1
            gpstruct.cov = {@covPERiso, gpstruct.cov};
        else
            gpstruct.cov = {@covPPERard_fast, gpstruct.cov};
        end
    else
        gpstruct.ncovoffset = 0;
    end

    % Prior and bounds on covariance length scale(s)
    covrange = (optimState.UB - optimState.LB)./optimState.scale;
    % Maximum length scale (normalized units, with PLB=-1 and PUB=1)
    covrange = min(100, 10*covrange);
    for i = 1:gpstruct.ncovlen
        % Wide prior on (log) normalized length scale 
        gpstruct.prior.cov{end+1} = {gaussPriorFunc, -1, 2^2};
        gpstruct.bounds.cov{end+1} = [log(TolMesh); log(covrange(i))];
    end
    
    % Prior and bounds on signal variance
    sf = exp(1);
    gpstruct.prior.cov{end+1} = {gaussPriorFunc, log(sf), 2^2};
    gpstruct.bounds.cov{end+1} = [log(TolFun); log(1e6*TolFun/TolMesh)];

    % Exponent of rational quadratic kernel
    if strcmpi(covtype,'rq')
        if numel(covextras) > 1; rqpriormean = covextras(2); else rqpriormean = 1; end
        gpstruct.prior.cov{end+1} = {gaussPriorFunc, rqpriormean, 1^2};
        gpstruct.bounds.cov{end+1} = [-5;5];
    end
    
    % Rotate GP axes (unsupported)
    if options.RotateGP && isfield(optimState,'C')
        gpstruct.C = inv(optimState.C);
    end
        
    %% GP likelihood function

    gpstruct.hyp.lik = [];    
    gpstruct.prior.lik = [];
    gpstruct.bounds.lik = [];
    
    % Warped likelihood (unsupported)
    if options.WarpFunc > 0
        error('Warped likelihoods not supported at the moment.');
        
        gpstruct.lik = {@likGaussWarpExact,{@warpPower}};
        % Add warped likelihood hyperparameters
        for iWarp = 1:options.WarpFunc
            gpstruct.hyp.lik(end+1) = 0;
            % gpstruct.prior.lik{end+1} = {gaussPriorFunc, 0, 100^2};
            gpstruct.prior.lik{end+1} = {@priorDelta};
            gpstruct.hyp.lik(end+1) = 0;
            gpstruct.prior.lik{end+1} = {gaussPriorFunc, log(1), 0.1^2};
            gpstruct.bounds.lik{end+1} = [-Inf; Inf];
            gpstruct.bounds.lik{end+1} = [-Inf; Inf];        
        end
    else
        % Default Gaussian likelihood
        gpstruct.lik = @likGauss;
    end
    
    if ~isempty(gplik)   % Known noise level
        error('Fixed noise not supported.');
        gpstruct.prior.lik{end+1} = {@delta_prior};
        gpstruct.hyp.lik(end+1) = gplik;
        gpstruct.knownNoise = true;

    else                % Unknown noise level (even for deterministic functions, helps regularization)
        % likmu = log(1); liks2 = 1^2;
        likmu = log(options.NoiseSize(1));
        if numel(options.NoiseSize) > 1 && isfinite(options.NoiseSize(2))
            liks2 = options.NoiseSize(2)^2;
        else
            liks2 = 1^2;
        end
        gpstruct.hyp.lik(end+1) = likmu;
        gpstruct.prior.lik{end+1} = {gaussPriorFunc, likmu, liks2};
        gpstruct.knownNoise = false;
    end
        
    gpstruct.hyp.lik = gpstruct.hyp.lik(:);

    % Bounds on likelihood noise parameter
    gpstruct.bounds.lik{end+1} = [log(TolFun)-1; 5];

    %% GP mean function
    gpstruct.hyp.mean = 0;

    gpstruct.mean = @meanConst;             % Constant mean function
    if options.gpFixedMean
        gpstruct.prior.mean = {{@priorDelta}};  % Fixed mean        
    else
        gpstruct.prior.mean = {{gaussPriorFunc,0,1^2}}; 
    end
    gpstruct.bounds.mean{1} = [-Inf; Inf];

    %% GP sampling weight
    Nsamples = max(1,options.gpSamples);
    gpstruct.hypweight = ones(1,Nsamples)/Nsamples;
    gpstruct.hypmean = [];

    %% Finalize GP
        
    % Verify every field
    gpstruct = gpset(gpstruct);
    
    % Initial length scale
    gpstruct.lenscale = 1;    
    gpstruct.pollscale = ones(1, D);
    
    % GP effective length scale radius
    gpstruct.effectiveradius = 1;    
        
    % Initial GP variability scale
    gpstruct.sf = sf;
    
    % Do not marginalize hyperparameters (needs gpml_extensions toolbox)
    gpstruct.marginalize = 0;
            
else    % Update existing GP struct
    
    % We assume that the mean of the GP is higher than what we see
    ymean = prctile1(gpstruct.y,options.gpMeanPercentile);
    yrange = feval(options.gpMeanRangeFun, ymean, gpstruct.y);
    
    %% Update empirical prior for GP likelihood
    
    % Likelihood prior (scales with mesh size)
    gpstruct.prior.lik{end}{2} = log(options.NoiseSize(1)) + options.MeshNoiseMultiplier*log(MeshSize);
    % gpstruct.prior.lik{end}{2} = min(log(TolFun),log(MeshSize));
    
    if options.WarpFunc > 0     % Warped likelihood (unsupported)
        ywarp = prctile1(gpstruct.y,50) + 5*(prctile1(gpstruct.y,50) - prctile1(gpstruct.y,15.87));        
        for i = 1:numel(gpstruct.hyp)
            gpstruct.hyp(i).lik(1) = ywarp;
        end
    end    
    
    %% Update empirical prior for GP mean
    
    gpstruct.prior.mean{1}{2} = ymean;
    if ~options.gpFixedMean
        gpstruct.prior.mean{1}{3} = yrange.^2/4;
    end
    % gpstruct.bounds.mean{1}(1) = min(gpstruct.y);
    
    for i = 1:length(gpstruct.hyp)
        if options.gpFixedMean
            gpstruct.hyp(i).mean = ymean;
        end
        % gpstruct.hyp(i).mean = max(gpstruct.hyp(i).mean, gpstruct.bounds.mean{1}(1));
        % gpstruct.hyp(i).mean = ymean;
    end    
    
    %% Update empirical prior for GP covariance
    
    % Compute max and min distance between training inputs
    temp(1,:,:) = gpstruct.x';
    
    switch lower(options.gpCovPrior)
        
        case 'iso'    
            dist = squeeze(udist(gpstruct.x,temp,1,optimState));
            dist(dist == 0) = NaN;
            uu = 0.5*log(max(dist(:)));
            ll = 0.5*log(min(dist(:)));

            % Empirical prior on covariance lengths
            covmu = 0.5*(uu+ll);
            covsigma = (uu-ll)/2;

            for i = (1:gpstruct.ncovlen) + gpstruct.ncovoffset
                gpstruct.prior.cov{i} = {gaussPriorFunc,covmu,covsigma^2};
            end
        
        case 'ard'
            dist = udist(gpstruct.x,temp,1,optimState,1);
            dist(dist == 0) = NaN;
            uu = 0.5*log(squeeze(max(max(dist,[],1),[],3)));
            ll = 0.5*log(squeeze(min(min(dist,[],1),[],3)));

            % Empirical prior on covariance lengths
            covmu = 0.5*(uu+ll);
            covsigma = (uu-ll)/2;
            % covmu =  0.5*log(nanmean(dist(:)));

            % [covmu covsigma 0.5*log(nanmean(dist(:)))]

            for i = (1:gpstruct.ncovlen) + gpstruct.ncovoffset
                 gpstruct.prior.cov{i} = {gaussPriorFunc, ...
                     0.5*(mean(covmu)+covmu(i)),0.5*(mean(covsigma.^2) + covsigma(i)^2)};
                 % gpstruct.prior.cov{i} = {gaussPriorFunc,covmu,covsigma^2};
            end
        
        otherwise
            error('Unknown prior for covariance parameters.');
    end    
    
    % Adjust prior length scales for periodic variables (mapped to unit circle)
    if any(optimState.periodicvars)                
        per = log((optimState.UB - optimState.LB)./optimState.scale);
        for d = find(optimState.periodicvars)
            gpstruct.prior.cov{d + gpstruct.ncovoffset}{2} = ...
                gpstruct.prior.cov{d + gpstruct.ncovoffset}{2} - per(d); % + log(2*pi);
        end
    end

    % Empirical prior on signal variance
    if options.WarpFunc > 0
        warp = gpstruct.lik{2};
        ng = feval(warp{:});
        gy = feval(warp{:},gpstruct.y,gpstruct.hyp.lik(1:ng));
        sdy = log(std(gy));
    else
        sdy = log(std(gpstruct.y));
    end
    gpstruct.prior.cov{gpstruct.ncovoffset+gpstruct.ncovlen+1} = {gaussPriorFunc,sdy,2^2};

    % Exponent of rational quadratic kernel
    %if strcmpi(covtype,'rq')
    %     gpstruct.prior.cov{gpstruct.ncovlen+2}{2} = log(optimState.funccount);
    %end    

    %plot(log(diff(sort(gpstruct.y))));
    %drawnow;
    
end

%% GP inference method

gpstruct.infMethod = 'exact';

if options.WarpFunc > 0
    gpstruct.inf = {@infPrior, @infExactWarp, gpstruct.prior};        
else
    gpstruct.inf = {@infPrior_fast, {@infExact_fastrobust,options.CholAttempts}, gpstruct.prior};
end