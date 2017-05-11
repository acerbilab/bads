function gpstruct = gpdefStationaryNew(covtype,covtheta,D,gplik,optimState,options,gpstruct,gpmlext)
%BGA_GPINIT

if nargin < 7; gpstruct = []; end
if nargin < 8 || isempty(gpmlext); gpmlext = 0; end     % Use gpml_extensions

MeshSize = optimState.meshsize;
TolFun = options.TolFun;
TolMesh = optimState.TolMesh;

if gpmlext
    gaussPriorFunc = @gaussian_prior;
else
    gaussPriorFunc = @priorGauss;
end

covard = covtheta(1);   % First covariance parameter (flag)

% Initialize new gp struct
if isempty(gpstruct)

    %% gp covariance function
    gpstruct.covtype = lower(covtype);

    if gpmlext  % Covariance functions with extended API for Hessian calculations
        
        switch lower(covtype)
            case 'se'
                if covard;  gpstruct.cov = {@ard_sqdexp_covariance};
                else        gpstruct.cov = {@isotropic_sqdexp_covariance}; end
            case 'matern5' % Not supported yet
                if covard;  gpstruct.cov = {@ard_matern_covariance,5};
                else        gpstruct.cov = {@isotropic_matern_covariance,5}; end
            case 'rq'
%                 if covard;  gpstruct.cov = {@ard_ratquad_covariance};
                if covard;  gpstruct.cov = {@ard_ratquad_covariance_fast};
                else        error('Rational quadratic covariance not supported yet in gpml extensions.');   end
            otherwise
                error(['Covariance type ''' covtype ''' unavailable in gpml_extensions.']);
        end
        
    else        % All stationary covariance functions
        
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
    end
    
    ncov = eval(feval(gpstruct.cov{:}));
    gpstruct.hyp.cov = zeros(ncov,1);
    
    if strcmpi(covtype,'rq'); ncovlen = ncov - 2; 
    else ncovlen = ncov - 1; end
    % gpstruct.ncov = ncov;
    gpstruct.ncovlen = ncovlen;
    
    gpstruct.prior.cov = [];
    gpstruct.bounds.cov = [];
    
    % Prior and bounds for covariance periodic length
    if any(optimState.periodicvars)
        per = log((optimState.UB - optimState.LB)./optimState.scale);
        
        % Infinite period for non-periodic dimensions
        per(~optimState.periodicvars) = Inf;
                
        gpstruct.hyp.cov = [per(:); gpstruct.hyp.cov];
        
        for i = 1:gpstruct.ncovlen
            gpstruct.prior.cov{end+1} = {@priorDelta};
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
    covrange = min(100, 10*covrange);
    for i = 1:gpstruct.ncovlen
        gpstruct.prior.cov{end+1} = {gaussPriorFunc, -1, 2^2};
        gpstruct.bounds.cov{end+1} = [log(TolMesh); log(covrange(i))]; % log(100)
    end
    
    % Prior and bounds on signal variance
    sf = exp(1);
    gpstruct.prior.cov{end+1} = {gaussPriorFunc, log(sf), 2^2};
    gpstruct.bounds.cov{end+1} = [log(TolFun); log(1e6*TolFun/TolMesh)];

    % Exponent of rational quadratic kernel
    if strcmpi(covtype,'rq')
        if numel(covtheta) > 1; rqpriormean = covtheta(2); else rqpriormean = 1; end
        gpstruct.prior.cov{end+1} = {gaussPriorFunc, rqpriormean, 1^2};
        gpstruct.bounds.cov{end+1} = [-5;5];
        % gpstruct.prior.cov{end+1} = {@priorDelta};
    end
    
    %% gp likelihood function

    if options.WarpFunc > 0
        gpstruct.lik = {@likGaussWarpExact,{@warpPower}};
    else
        gpstruct.lik = @likGauss;   % Gaussian likelihood
    end
    
    gpstruct.hyp.lik = [];    
    gpstruct.prior.lik = [];
    
    for iWarp = 1:options.WarpFunc
        gpstruct.hyp.lik(end+1) = 0;
        % gpstruct.prior.lik{end+1} = {gaussPriorFunc, 0, 100^2};
        gpstruct.prior.lik{end+1} = {@priorDelta};
        gpstruct.hyp.lik(end+1) = 0;
        gpstruct.prior.lik{end+1} = {gaussPriorFunc, log(1), 0.1^2};
    end
    
    if ~isempty(gplik)   % Known noise level
        error('a');
        gpstruct.prior.lik{end+1} = {@delta_prior};
        gpstruct.hyp.lik(end+1) = gplik;
        gpstruct.knownNoise = true;

    else                % Unknown noise level
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
    gpstruct.bounds.lik = [];
    if options.WarpFunc > 0
        gpstruct.bounds.lik{end+1} = [-Inf; Inf];
        gpstruct.bounds.lik{end+1} = [-Inf; Inf];        
    end
    gpstruct.bounds.lik{end+1} = [log(TolFun)-1; 5];

    %% gp mean function

    gpstruct.hyp.mean = 0;
    if gpmlext
        gpstruct.mean = @constant_mean;
        if options.gpFixedMean
            gpstruct.prior.mean = {{@delta_prior, 0}};  % Fixed mean        
        else
            gpstruct.prior.mean = {{gaussPriorFunc,0,1^2}}; 
        end
    else        
        gpstruct.mean = @meanConst;             % Constant mean function
        if options.gpFixedMean
            gpstruct.prior.mean = {{@priorDelta}};  % Fixed mean        
        else
            gpstruct.prior.mean = {{gaussPriorFunc,0,1^2}}; 
        end
    end
    gpstruct.bounds.mean{1} = [-Inf; Inf];

    %% gp sampling weight

    Nsamples = max(1,options.gpSamples);
    gpstruct.hypweight = ones(1,Nsamples)/Nsamples;
    gpstruct.hypmean = [];

    % Check every field
    gpstruct = gpset(gpstruct);
    
    % Initial length scale
    gpstruct.lenscale = 1;    
    gpstruct.pollscale = ones(1, D);
    
    % gp effective length scale radius
    gpstruct.effectiveradius = 1;    
        
    % Initial gp variability scale
    gpstruct.sf = sf;
    
    % Store whether using gpml_extensions toolbox
    gpstruct.marginalize = gpmlext;
        
else    % Update existing gp struct
    
    ymean = prctile1(gpstruct.y,options.gpMeanPercentile); % Write your own function here
    % yrange = prctile1(gpstruct.y,75) - prctile1(gpstruct.y,25);
    % yrange = (ymean - prctile1(gpstruct.y,50))/5*2;
    yrange = feval(options.gpMeanRangeFun, ymean, gpstruct.y);
    
    %% Update likelihood
    
    % Likelihood prior
    gpstruct.prior.lik{end}{2} = log(options.NoiseSize(1)) + options.MeshNoiseMultiplier*log(MeshSize);
    % gpstruct.prior.lik{end}{2} = min(log(TolFun),log(MeshSize));
    % warning('Ho modificato il prior su likelihood in gpdefStationary!');
    
    if options.WarpFunc > 0
        ywarp = prctile1(gpstruct.y,50) + 5*(prctile1(gpstruct.y,50) - prctile1(gpstruct.y,15.87));
        % gpstruct.prior.lik{1}{2} = ymean;
        % yrange = prctile1(gpstruct.y,90) - prctile1(gpstruct.y,10);
        % yrange = max(gpstruct.y) - min(gpstruct.y);
        % gpstruct.prior.lik{1}{3} = yrange.^2/4;
        %gpstruct.prior.lik{2}{2} = 0; % log(yrange.^2/4);
        %gpstruct.prior.lik{2}{3} = 0.5^2;
        
        for i = 1:numel(gpstruct.hyp)
            gpstruct.hyp(i).lik(1) = ywarp;
            % gpstruct.hyp(i).lik(1) = gpstruct.prior.lik{1}{2};        
            % gpstruct.hyp(i).lik(2) = gpstruct.prior.lik{2}{2};
        end
    end    
    
    %% Update mean
    
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
    
    %% Update covariance
    
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

    % Prior on signal variance
    if options.WarpFunc > 0
        warp = gpstruct.lik{2};
        ng = feval(warp{:});
        gy = feval(warp{:},gpstruct.y,gpstruct.hyp.lik(1:ng));
        sdy = log(std(gy));
    else
        sdy = log(std(gpstruct.y));
    end
    %sdy = 0;
    gpstruct.prior.cov{gpstruct.ncovoffset+gpstruct.ncovlen+1} = {gaussPriorFunc,sdy,2^2};

    % Exponent of rational quadratic kernel
    %if strcmpi(covtype,'rq')
    %     gpstruct.prior.cov{gpstruct.ncovlen+2}{2} = log(optimState.funccount);
    %end    

    %plot(log(diff(sort(gpstruct.y))));
    %drawnow;
    
end


%% gp inference method

% gpstruct.hyp.lik(:)'

gpstruct.infMethod = 'exact';

if gpmlext
    prior = @(theta_) independent_prior(gpstruct.prior,theta_);
    gpstruct.inf = {@inference_with_prior, @exact_inference_robust, prior};
    % gpstruct.inf = {@inference_with_prior, @exact_inference_fast, prior};
else
    if options.WarpFunc == 0
        gpstruct.inf = {@infPrior_fast, {@infExact_fastrobust,options.CholAttempts}, gpstruct.prior};
    else
        gpstruct.inf = {@infPrior, @infExactWarp, gpstruct.prior};        
    end
end

end
