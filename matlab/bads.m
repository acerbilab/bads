function [x,fval,exitflag,output,funValues,gpstruct] = bads(fun,x0,LB,UB,PLB,PUB,options,varargin)
%BADS Constrained optimization using Bayesian pattern search
%   BADS attempts to solve problems of the form:
%       min F(X)  subject to:  LB <= X <= UB
%        X
%
%   X = BADS(FUN,X0,LB,UB) starts at X0 and finds a local minimum X to
%   the function FUN. FUN accepts input X and returns a scalar function
%   value evaluated at X. X0 may be a scalar, vector or empty matrix. If X0 
%   is empty, the starting point is chosen from the initial mesh only.
%   LB and UB define a set of lower and upper bounds on the design 
%   variables, X, so that a solution is found in the range LB <= X <= UB.
%   LB and UB can be scalars or vectors. If scalars, the bound is 
%   replicated in each dimension. Use empty matrices for LB and UB if no 
%   bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; set 
%   UB(i) = Inf if X(i) is unbounded above. Note that: 
%      - if LB and/or UB contain unbounded variables, the respective 
%        values of PLB and/or PUB need to be specified (see below).
%      - if X0 is empty, PLB and PUB need to be specified as vectors.
%
%   X = BADS(FUN,X0,LB,UB,PLB,PUB) specifies a set of plausible lower and
%   upper bounds such that LB <= PLB <= X0 <= PUB <= UB. Both PLB and PUB
%   need to be finite. PLB and PUB are used to design the initial mesh of 
%   the pattern search, and represent a plausible range for the 
%   optimization variables. As a rule of thumb, set PLB and PUB such that 
%   there is ~ 95% probability that the minimum is found within the box 
%   (where in doubt, just set PLB=LB and PUB=UB).
%
%   X = BADS(FUN,X0,LB,UB,PLB,PUB,options) minimizes with the default 
%   optimization parameters replaced by values in the structure OPTIONS.
%   BPS('defaults') returns the default OPTIONS struct.
%
%   [X,FVAL] = BADS(...) returns FVAL, the value of the objective function 
%   FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = BADS(...) returns EXITFLAG which describes the exit 
%   condition of BADS. Possible values of EXITFLAG and the corresponding 
%   exit conditions are
%
%     0  Maximum number of function evaluations or iterations reached.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = BADS(...) returns a structure OUTPUT with the
%   following information:
%   [...] to be written [...]
%
%   [...]
%
%   Author: Luigi Acerbi
%   Release date: Mar 14, 2017
%   Version: 0.7-alpha -- INTERNAL USE ONLY, DO NOT DISTRIBUTE


%
% LAST CHANGES: 
% gpTrainingSet:         - uCheck set to projection = 1 line 334
%
% TO DO:
% - check options.FunValues provided as argument
% - fix restarts e multibayes
% - compute func running time and do more stuff if func is slow
% - understand if the warping is working correctly, test moar
% - testa se la search Hedge è ancora necessaria?
% - aggiungi warning se la mesh cerca di espandere oltre MaxPollGridNumber
%   (sintomo di misspecification dei bound)
% - in private, functions che iniziano con init si confondono con initSorbole


%% Default options

defopts.Display                 = 'iter                 % Level of display ("iter", "notify", "final", or "off")';
defopts.Plot                    = 'off                  % Show optimization plots ("profile", "scatter", or "off")';
defopts.Debug                   = 'off                  % Debug mode, plot additional info';

% Termination conditions
defopts.MaxIter                 = '200*nvars            % Max number of iterations';
defopts.MaxFunEvals             = '500*nvars            % Max number of objective fcn evaluations';
defopts.TolMesh                 = '1e-6                 % Tolerance on mesh size';
defopts.TolFun                  = '1e-3                 % Min significant change of objective fcn';
defopts.TolStallIters           = '4 + floor(nvars/2)   % Max iterations with no significant change (doubled under uncertainty)';

defopts.Ninit                   = 'nvars                % Number of initial objective fcn evaluations';
defopts.InitFcn                 = '@initSobol           % Initialization function';
% defoptions.InitFcn            = '@initLHS';
defopts.PollMethod              = '@pollMADS2N          % Poll function';
defopts.Nbasis                  = '200*nvars';

defopts.Restarts                = '0                    % Number of restart attempts';
defopts.CacheSize               = '1e4                  % Size of cache';
defopts.FunValues               = '[]                   % Struct with pregress fcn evaluations';
defopts.PeriodicVars            = '[]                   % Array with indices of periodic variables';

defopts.TolImprovement          = '1                    % Minimum significant improvement at unit mesh size';
defopts.ForcingExponent         = '3/2                  % Exponent of forcing function';
defopts.PollMeshMultiplier      = '2                    % Mesh multiplicative factor between iterations';
defopts.IncumbentSigmaMultiplier = '0.1                 % Multiplier to incumbent uncertainty for acquisition functions';
defopts.AlternativeIncumbent    = 'off                  % Use alternative incumbent offset';
defopts.AdaptiveIncumbentShift  = 'off                  % Adaptive multiplier to incumbent uncertainty';
defopts.NonlinearScaling        = 'on                   % Allow nonlinear rescaling of variables if deemed useful';
defopts.gpRescalePoll           = '1                    % GP-based geometric scaling factor of poll vectors';
defopts.FitnessShaping          = 'off                  % Nonlinear rescaling of objective fcn';
defopts.WarpFunc                = '0                    % GP warping function type';
defopts.OptimToolbox            = '[]                   % Use MATLAB Optimization Toolbox (if empty, determine at runtime)';

% Search properties
defopts.Nsearch                 = '2^12                 % Number of candidate search points';
defopts.Nsearchiter             = '2                    % Number of optimization iterations for search';
defopts.ESbeta                  = '1                    % Multiplier in ES';
defopts.ESstart                 = '0.25                 % Starting scale value in ES';
defopts.SearchImproveFrac       = '0                    % Fraction of candidate search points with (slower) improved estimate';
defopts.SearchScaleSuccess      = 'sqrt(2)              % Search radius expansion factor for successful search';
defopts.SearchScaleIncremental  = '2                    % Search radius expansion factor for incremental search';
defopts.SearchScaleFailure      = 'sqrt(0.5)            % Search radius contraction factor for failed search';
defopts.SearchFactorMin         = '0.5';
% defopts.SearchMethod            = '{@searchWCMscale,1}  % Search function(s)';
defopts.SearchMethod            = '{@searchHedge,{{@searchES,1,1},{@searchES,2,1}}}  % Search function(s)';
defopts.SearchGridNumber        = '10                   % iteration scale factor between poll and search';
defopts.MaxPollGridNumber       = '0                    % Maximum poll integer';
defopts.SearchGridMultiplier    = '2                    % multiplier integer scale factor between poll and search';
defopts.SearchSizeLocked        = 'on                   % Relative search scale factor locked to poll scale factor';
% defopts.SearchNtry              = 'max(2*nvars,5+nvars) % Number of searches per iteration';
defopts.SearchNtry              = 'max(nvars,floor(3+nvars/2)) % Number of searches per iteration';
defopts.SearchMeshExpand        = '0                    % Search-triggered mesh expansion after this number of successful search rounds';
defopts.SearchMeshIncrement     = '1                    % Mesh size increment after search-triggered mesh expansion';
defopts.SearchOptimize          = 'no                   % Further optimize acquisition function';

% Gaussian Process properties
defopts.Ndata                   = '50 + 10*nvars       % Number of training data (doubled under uncertainty)';
% defopts.Ndata                   = 'min(max(50,10*nvars),200)    % Number of training data (doubled under uncertainty)';
defopts.MinNdata                = '50                   % Minimum number of training data (doubled under uncertainty)';
defopts.BufferNdata             = '100                  % Max number of training data removed if too far from current point';
defopts.gpSamples               = '0                    % Hyperparameters samples (0 = optimize)';
defopts.gpMarginalize           = 'off                  % Approximate marginalization of hyper-parameters';
defopts.MinRefitTime            = '2*nvars              % Minimum fcn evals before refitting the gp';
defopts.PollTraining            = 'yes                  % Train gp also during poll step';
defopts.DoubleRefit             = 'off                  % Try a second fit';
defopts.gpMeanPercentile        = '90                   % Percentile of empirical GP mean';
defopts.gpMeanRangeFun          = '@(ym,y) (ym - prctile(y,50))/5*2   % Empirical range of hyperprior over the mean';
defopts.gpdefFcn                = '{@gpdefStationaryNew,''rq'',[1,1]}  % GP definition fcn';
%defopts.gpdefFcn                = '{@gpdefStationaryNew,''matern5'',1}     % GP definition fcn';
defopts.gpMethod                = 'grid                 % GP fit method';
defopts.gpCluster               = 'no                   % Cluster additional points during training';
defopts.RotateGP                = 'no                   % Rotate GP basis';
defopts.gpRadius                = '3                    % Radius of training set';
defopts.UseEffectiveRadius      = 'yes                   %';
defopts.gpCovPrior              = 'iso                  % GP hyper-prior over covariance';
defopts.gpFixedMean             = 'no';
defopts.FitLik                  = 'yes                  % Fit the likelihood term';
defopts.PollAcqFcn              = '@acqNegEI            % Acquisition fcn for poll step';
defopts.SearchAcqFcn            = '@acqNegEI            % Acquisition fcn for search step';
defopts.AcqHedge                = 'off                  % Hedge acquisition function';
defopts.CholAttempts            = '0                    % Attempts at performing the Cholesky decomposition';
defopts.NoiseNudge              = '[1 0]                % Increase nudge to noise in case of Cholesky failure';
defopts.RemovePointsAfterTries  = '1                    % Start removing training points after this number of failures';

defopts.UncertaintyHandling     = 'off                  % Explicit handling of noise';
defopts.NoiseObj                = 'off                  % Objective function returns estimated noise value as second argument';
defopts.UncertainIncumbent      = 'yes                  % Treat incumbent as if uncertain regardless of uncertainty handling';
defopts.NoiseSize               = 'sqrt(options.TolFun) % Base observation noise magnitude';
defopts.MeshNoiseMultiplier     = '0.5                  % Contribution to log noise magnitude from log mesh size (0 for noisy functions)';
defopts.TolPoI                  = '1e-6/nvars           % Threshold probability of improvement (PoI); set to 0 to always complete polling';
defopts.SkipPoll                = 'yes                  % Skip polling if PoI below threshold, even with no success';
defopts.ConsecutiveSkipping     = 'yes                  % Allow consecutive incomplete polls';
defopts.SkipPollAfterSearch     = 'yes                  % Skip polling after successful search';
% defopts.MinFailedPollSteps      = 'ceil(sqrt(nvars))    % Number of failed fcn evaluations before skipping is allowed';
defopts.MinFailedPollSteps      = 'Inf                  % Number of failed fcn evaluations before skipping is allowed';
defopts.NormAlphaLevel          = '1e-6                 % Alpha level for normality test of gp predictions';
defopts.AccelerateMesh          = 'on                   % Accelerate mesh contraction';
defopts.AccelerateMeshSteps     = '3                    % Accelerate mesh after this number of stalled iterations';
defopts.SloppyImprovement       = 'yes                  % Move incumbent even after insufficient improvement';
defopts.HessianUpdate           = 'no                   % Update Hessian as you go';
defopts.HessianMethod           = 'bfgs                 % Hessian update method';
defopts.HessianAlternate        = 'no                   % Alternate Hessian iterations';

defopts.HedgeGamma              = '0.125';
defopts.HedgeBeta               = '1e-3/options.TolFun';
defopts.HedgeDecay              = '0.1^(1/(2*nvars))';

defopts.TrueMinX                = '[]                   % Location of the global minimum (for visualization only)';

%% If called with no arguments or with 'defaults', return default options
if nargin < 1 || strcmpi(fun, 'defaults')
    if nargin < 1
        fprintf('Default options returned (type "help bads" for help).\n');
    end
    x = defopts;
    return;
end


%% Check that all required subfolders are on the MATLAB path

subfolders = {'acq','gpdef','gpml_fast','init','poll','search','utils','warp'};
pathCell = regexp(path, pathsep, 'split');
baseFolder = fileparts(mfilename('fullpath'));

onPath = true;
for iFolder = 1:numel(subfolders)
    folder = [baseFolder,filesep,subfolders{iFolder}];    
    if ispc  % Windows is not case-sensitive
      onPath = onPath & any(strcmpi(folder, pathCell));
    else
      onPath = onPath & any(strcmp(folder, pathCell));
    end
end

% ADDPATH is slow, call it only if folders are not on path
if ~onPath
    addpath(genpath(fileparts(mfilename('fullpath'))));
end

%% Initialize display printing options

if ~isfield(options,'Display') || isempty(options.Display)
    options.Display = defopts.Display;
end

switch lower(options.Display(1:3))
    case {'not','notify','notify-detailed'}
        prnt = 1;
    case {'non','none','off'}
        prnt = 0;
    case {'ite','all','iter','iter-detailed'}
        prnt = 3;
    case {'fin','final','final-detailed'}
        prnt = 2;
    otherwise
        prnt = 1;
end

%% Initialize variables and algorithm structures

if nargin < 3 || isempty(LB); LB = -Inf; end
if nargin < 4 || isempty(UB); UB = Inf; end
if nargin < 5; PLB = []; end
if nargin < 6; PUB = []; end
if nargin < 7; options = []; end

if isempty(x0)
    if prnt > 2
        fprintf('X0 not specified. Taking the number of dimensions from PLB and PUB...');
    end
    if isempty(PLB) || isempty(PUB)
        error('If no starting point is provided, PLB and PUB need to be specified.');
    end    
    x0 = NaN(size(PLB));
    if prnt > 2
        fprintf(' NVARS = %d.\n', numel(x0));
    end
    SkipInitPoint = 1;  % No starting point provided
else
    SkipInitPoint = 0;  % Starting point provided
end

nvars = numel(x0);
optimState = [];

% Initialize algorithm options
options = initoptions(nvars,defopts,options);    

% Initalize and transform variables
[u0,LB,UB,PLB,PUB,MeshSizeInteger,MeshSize,TolMesh,optimState] = ...
    initvars(x0,LB,UB,PLB,PUB,optimState,options);
options.TolMesh = TolMesh;
    
optimState = updateSearchBounds(optimState);
optimState.searchfactor = 1;
optimState.sdlevel = options.IncumbentSigmaMultiplier;
optimState.searchcount = options.SearchNtry;    % Skip search at first iteration

% Create vector of ES weights (only for searchES)
es_iter = options.Nsearchiter;
es_mu = options.Nsearch/es_iter;
es_lambda = es_mu;
optimState.es = ESupdate(es_mu,es_lambda,es_iter);

% Store objective function
if ischar(fun); fun = str2func(fun); end
optimState.fun = fun;
if isempty(varargin)
    funwrapper = fun;   % No additional function arguments passed
else
    funwrapper = @(u_) fun(u_,varargin{:});
end

% Initialize function logger
[~,optimState] = funlogger([],u0,optimState,'init',options.CacheSize,options.NoiseObj);

%% Initial function evaluations

iter = 0;
optimState.iter = iter;

if prnt > 2
    if options.UncertaintyHandling
        displayFormat = ' %5.0f       %5.0f    %12.6g    %12.6g    %12.6g    %20s    %s\n';
        fprintf(' Iteration    f-count      E[f(x)]        SD[f(x)]      MeshScale          Method          Actions\n');
    else
        displayFormat = ' %5.0f       %5.0f    %12.6g    %12.6g    %20s    %s\n';
        fprintf(' Iteration    f-count          f(x)          MeshScale          Method          Actions\n');
    end
else
    displayFormat = [];
end

% Evaluate starting point and initial mesh
[u,fval,isFinished,optimState] = ...
    initmesh(u0,funwrapper,SkipInitPoint,optimState,options,prnt,displayFormat);
    
if options.UncertaintyHandling  % Current uncertainty in estimate
    fsd = options.NoiseSize(1);
else
    fsd = 0;
end
optimState.fsd = fsd;

ubest = u;                      % Current best minumum location
optimState.usuccess = ubest;    % Store sequence of successful x and y values
optimState.fsuccess = fval;
optimState.u = u;

% Initialize Gaussian Process (GP) structure
if options.FitLik; gplik = []; else gplik = log(options.TolFun); end
gpstruct = feval(options.gpdefFcn{:},nvars,gplik,optimState,options,[],options.gpMarginalize);
gpstruct.fun = funwrapper;
fhyp = gpstruct.hyp;
if options.RotateGP && isfield(optimState,'C'); gpstruct.C = inv(optimState.C); end

fvalhistory = Inf(min(1000,options.MaxIter),1);

% Initialize struct with GP prediction statistics
optimState.gpstats = savegpstats([],[],[],[],ones(1,max(1,options.gpSamples)));

% List of points at the end of each iteration
optimState.iterList.u = [];
optimState.iterList.fval = [];
optimState.iterList.fsd = [];
optimState.iterList.fhyp = [];

optimState.lastfitgp = -Inf;    % Last fcn evaluation for which the gp was trained
lastskipped = 0;                % Last skipped iteration
pollCount = 0;
optimState.hedge = [];

SearchSuccesses = 0;
SearchSpree = 0;
Restarts = options.Restarts;

%% Optimization loop
iter = 1;

while ~isFinished
    optimState.iter = iter;
    refitted = 0;                       % gp refitted this iteration
    action = [];
                
    MeshSize = options.PollMeshMultiplier^MeshSizeInteger;
    if options.SearchSizeLocked
        optimState.SearchSizeInteger = min(0,MeshSizeInteger*options.SearchGridMultiplier - options.SearchGridNumber);
    end
    optimState.meshsize = MeshSize;
    optimState.searchmeshsize = options.PollMeshMultiplier.^optimState.SearchSizeInteger;
    
    % Update bounds for search mesh
    optimState = updateSearchBounds(optimState);
    
    % Minimum improvement for a poll/search to be considered successful
    SufficientImprovement = options.TolImprovement*(MeshSize^options.ForcingExponent);
    if options.SloppyImprovement
        SufficientImprovement = max(SufficientImprovement, options.TolFun);
    end

    SearchSufficientImprovement = SufficientImprovement;
    % Multiple successful searches raise the bar for improvement
    % SearchSufficientImprovement = SufficientImprovement*(2^SearchSuccesses);
    optimState.SearchSufficientImprovement = SearchSufficientImprovement;
        
    %----------------------------------------------------------------------
    %% Search step

    tic
            
    if optimState.searchcount < options.SearchNtry && size(gpstruct.y,1) > 1
                
        % Check whether it is time to refit the GP
        [refitgp,UnrelGP,optimState] = IsRefitTime(optimState,options);
        if refitgp; gpstruct.post = []; end
        
        if isempty(gpstruct.post) % || optimState.funccount > options.Ndata
            % Local GP approximation on current point
            gpstruct = gpTrainingSet(gpstruct, ...
                options.gpMethod, ...
                u, ...
                [], ...    %             [upoll; gridunits(x,optimState)], ...
                optimState, ...
                options, ...
                refitgp);
        end
        
        % Update gp prediction at incumbent
        optimState = UpdateIncumbentPrediction(ubest,fhyp,optimState,gpstruct,options);        
        
        % Generate search set (normalized coordinates)        
        %tic
        optimState.searchcount = optimState.searchcount + 1;        
        [usearchset,optimState] = feval(options.SearchMethod{:}, ...
            u, ...
            gpstruct, ...
            LB, ...
            UB, ...
            optimState, ...
            options);    
        %toc
        
        % Enforce periodicity
        usearchset = periodCheck(usearchset,LB,UB,optimState);
        
        % Force candidate points on search grid
        usearchset = force2grid(usearchset, optimState);
                
        % Remove already evaluated or unfeasible points from search set
        usearchset = uCheck(usearchset,options.TolMesh,optimState,1);
                                
        if ~isempty(usearchset)

            ns = size(usearchset, 1);            
            ymu = zeros(numel(gpstruct.hyp),ns);
            ys = zeros(numel(gpstruct.hyp),ns);            

            % optimState = UpdateIncumbentPrediction(ubest,fhyp,optimState,gpstruct,options);
            
            try                
                if options.AcqHedge
                    [optimState.hedge,acqIndex,ymu,ys] = ...
                        acqPortfolio('acq',optimState.hedge,usearchset,optimState.ftarget,fstarget,gpstruct,optimState,options,SufficientImprovement);
                    index = acqIndex(optimState.hedge.chosen);
                    z = 1;
                else
                    % Batch evaluation of acquisition function on search set
                    [z,~,ymu,ys] = ...
                        feval(options.SearchAcqFcn{:},usearchset,optimState.ftarget,gpstruct,optimState,0);

                    % Evaluate best candidate point in original coordinates
                    [~,index] = min(z);
                end
            catch
                % Failed evaluation of the acquisition function
                index = []; z = [];
            end
                        
            % Randomly choose index if something went wrong
            if isempty(index) || ~isfinite(index); index = randi(size(usearchset,1)); end
            
            % Do something else if expected improvement is not big enough?
            if 0 && abs(z(index)) < SufficientImprovement
                abs(z(index))
                
                optimState2 = optimState;
                options2 = options;
                optimState2.searchfactor = 10;
                options2.SearchAcqFcn = {@acqLCB,[]};               
                usearchset = feval(@searchES,1,1, ...
                    u, ...
                    gpstruct, ...
                    LB, ...
                    UB, ...
                    optimState2, ...
                    options2);

                % Enforce periodicity
                usearchset = periodCheck(usearchset,LB,UB,optimState);

                % Force candidate points on search grid
                usearchset = force2grid(usearchset, optimState);

                % Remove already evaluated or unfeasible points from search set
                usearchset = uCheck(usearchset,options.TolMesh,optimState,1);
                
                [~,index] = min(z);
            end
            
            acqu = [];
            
            % Local optimization of the acquisition function 
            % (Note that generally it does not improve results)
            if options.SearchOptimize
                acqoptoptions = optimset('Display','off','GradObj','off','DerivativeCheck','off',...
                    'TolX',options.TolMesh,'TolFun',options.TolFun);

                try
                    acqu = usearchset(index,:);
                    acqoptoptions.MaxFunEval = options.Nsearch;
                    acqoptf = @(u_) feval(options.SearchAcqFcn{:},u_,optimState.ftarget,gpstruct,optimState,0);
                    % Limit seach within search box
                    % NEEDS TO BE ADJUSTED FOR PERIODIC VARIABLES
                    acqlb = max(LB, min(acqu,u - optimState.searchfactor*MeshSize));
                    acqub = min(UB, max(acqu,u + optimState.searchfactor*MeshSize));
                    [acqu,facq,~,output] = fmincon(acqoptf, ...
                        acqu,[],[],[],[],acqlb,acqub,[],acqoptoptions);
                    acqu = force2grid(acqu, optimState);
                    acqu = periodCheck(acqu,LB,UB,optimState,1);
                    acqu = uCheck(acqu,options.TolMesh,optimState,1);
                catch
                    acqu = [];
                end
            end
            
            if isempty(acqu); acqu = usearchset(index,:); end
            
            usearch = acqu;
            
            % Evaluate function on search point
            [fsearch,optimState] = funlogger(funwrapper,usearch,optimState,'iter');
            
            if ~isempty(z)
                % Save statistics of gp prediction
                optimState.gpstats = ...
                    savegpstats(optimState.gpstats,fsearch,ymu(:,index),ys(:,index),gpstruct.hypweight);
            end
            
            % Add search point to training set
            if ~isempty(usearch) && optimState.searchcount < options.SearchNtry
                gpstruct = gpTrainingSet(gpstruct, ...
                    'add', ...
                    usearch, ...
                    fsearch, ...
                    optimState, ...
                    options, ...
                    0);
            end
            
            if options.UncertaintyHandling
                gpstructnew = gpTrainingSet(gpstruct, ...
                    options.gpMethod, ...
                    usearch, ...
                    [], ...
                    optimState, ...
                    options, ...
                    0);
                
                % Compute estimated function value at point
                [~,~,fsearch,fs2] = gppred(usearch,gpstructnew);
                fsearchsd = sqrt(fs2);
            else
                fsearchsd = 0;
            end
            
            % Compute distance of search point from current point
            searchdist = sqrt(udist(ubest,usearch,gpstruct.lenscale,optimState));            
        else
            fsearch = fval;
            fsearchsd = 0;
            searchdist = 0;
        end
        
        % CMA-ES like estimation of local covariance structure
        if options.HessianUpdate && strcmpi(options.HessianMethod,'cmaes')
            optimState = covmatadapt(u,LB,UB,gpstruct,optimState,options);
        end
                
        % Evaluate search
        SearchImprovement = fval - fsearch + 0*(fsd  - fsearchsd);
        fvalold = fval;
        
        if (SearchImprovement > 0 && options.SloppyImprovement) ...
                || SearchImprovement > SearchSufficientImprovement
            if options.AcqHedge
                method = optimState.hedge.str{optimState.hedge.chosen};
            else
                method = feval(options.SearchMethod{:},[],[],[],[],optimState);
            end
            if SearchImprovement > SearchSufficientImprovement
                SearchSuccesses = SearchSuccesses + 1;
                searchstring = ['Successful search (' method ')'];
                optimState.usuccess = [optimState.usuccess; usearch];
                optimState.fsuccess = [optimState.fsuccess; fsearch];
                searchstatus = 'success';                
            else
                searchstring = ['Incremental search (' method ')'];
                % searchstring = [];
                searchstatus = 'incremental';
            end
            
            % Update incumbent point
            [ubest,fval,fsd,optimState,gpstruct] = UpdateIncumbent(ubest,fval,fsd,usearch,fsearch,fsearchsd,optimState,gpstruct,options);
            
            if options.UncertaintyHandling; gpstruct = gpstructnew; end
            gpstruct.post = []; % Reset posterior
        else
            searchstatus = 'failure';
            searchstring = [];
            
            % Decay of local curvature estimate
            %if options.HessianUpdate
            %    c1 = 0.1/nvars;
            %    rescaledLenscale = gpstruct.pollscale;
            %    rescaledLenscale = rescaledLenscale/sqrt(sum(rescaledLenscale.^2));                
            %    optimState.Binv = optimState.Binv*(1-c1) + c1*diag(rescaledLenscale.^2);
            %end            
        end
                        
        % Update portfolio acquisition function
        if options.AcqHedge && ~isempty(usearchset)            
            optimState.hedge = ...
                acqPortfolio('update',optimState.hedge,usearchset(acqIndex,:),fsearch,fsearchsd,gpstruct,optimState,options,SufficientImprovement,fvalold,MeshSize);            
        end
        
        % Update search portfolio (needs improvement)
        if ~options.AcqHedge && ~isempty(optimState.hedge) && ~isempty(usearch)
            optimState.hedge = ...
                acqPortfolio('update',optimState.hedge,usearch,fsearch,fsearchsd,gpstruct,optimState,options,SufficientImprovement,fvalold,MeshSize);            
        end
                
        % Update search statistics and search scale factor
       optimState = UpdateSearch(optimState,searchstatus,searchdist,options);
                
        if ~isempty(searchstring) && any(strcmpi(options.Display,{'all','iter'}))
            if options.UncertaintyHandling
                fprintf(displayFormat,iter,optimState.funccount,fval,fsd,MeshSize,searchstring,'');                
            else
                fprintf(displayFormat,iter,optimState.funccount,fval,MeshSize,searchstring,'');
            end
        end
            
    end % Search step

    % Decide whether to perform the poll step
    switch optimState.searchcount
        case {0, options.SearchNtry}    % Skipped or just finished search             
            optimState.searchcount = 0;
            if SearchSuccesses > 0 && options.SkipPollAfterSearch
                DoPollStep = 0;
                SearchSpree = SearchSpree + 1;
                if options.SearchMeshExpand > 0 && ...
                        mod(SearchSpree,options.SearchMeshExpand) == 0
                    MeshSizeInteger = min(MeshSizeInteger + options.SearchMeshIncrement, options.MaxPollGridNumber);
                end
            else
                DoPollStep = 1;
                SearchSpree = 0;
            end
            SearchSuccesses = 0;
            % optimState.searchfactor = 1;
        otherwise                       % In-between searches, no poll
            DoPollStep = 0;
    end
            
    %----------------------------------------------------------------------
    %% Poll step  

    u = ubest;
    
    if DoPollStep
        
        PollImprovement = 0;            % Improvement so far
        upollbest = u;                  % Best poll point
        fpollbest = fval;               % Objective value at best point
        fpollhyp = fhyp;                % gp hyper-parameters at best point
        fpollbestsd = fsd;              % Uncertainty of objective func
        optimState.pollcount = 0;       % Poll iterations
        goodpoll = 0;                   % Found a good poll
        B = [];                         % Poll basis
        upoll = [];                     % Poll vectors
        unew = [];

        % Poll loop
        while (~isempty(upoll) || isempty(B)) ...
                && optimState.funccount < options.MaxFunEvals

            % Fill in basis vectors
            Bnew = feval(options.PollMethod{:}, ...
                B, ...
                u, ...
                gpstruct, ...
                LB, ...
                UB, ...
                optimState, ...
                options);

            % Create new poll vectors
            if ~isempty(Bnew)
                % GP-based vector scaling
                vv = bsxfun(@times,Bnew*MeshSize,gpstruct.pollscale);
                
                upollnew = periodCheck(bsxfun(@plus,u,vv),LB,UB,optimState);
                upollnew = uCheck(upollnew,options.TolMesh,optimState,0);
                upoll = [upoll; upollnew];
                B = [B; Bnew];                    
            end
            
            % Cannot refill poll vector set, stop polling
            if isempty(upoll); break; end
            
            % Check whether it is time to refit the GP
            [refitgp,UnrelGP,optimState] = IsRefitTime(optimState,options);
            if ~options.PollTraining && iter > 1; refitgp = 0; end
            if refitgp || optimState.pollcount == 0; gpstruct.post = []; end
            
            % Local GP approximation around polled points
            if isempty(gpstruct.post) % || optimState.funccount > options.Ndata
                gpstruct = gpTrainingSet(gpstruct, ...
                    options.gpMethod, ...
                    u, ...
                    upoll, ...
                    optimState, ...
                    options, ...
                    refitgp);
                if refitgp; refitted = 1; end
            end

            gpstruct.ystar = fpollbest;     % Best reference value

            optimState = UpdateIncumbentPrediction(upollbest,fpollhyp,optimState,gpstruct,options);
                        
            % Evaluate acquisition function on poll vectors
            if options.AcqHedge
                [optimState.hedge,acqIndex,ymu,ys,fm,fs] = ...
                    acqPortfolio('acq',optimState.hedge,upoll,optimState.ftarget,fstarget,gpstruct,optimState,options,SufficientImprovement);
                index = acqIndex(optimState.hedge.chosen);
            else
                % Batch evaluation of acquisition function on search set
                [z,~,ymu,ys,fm,fs] = feval(options.PollAcqFcn{:},upoll,optimState.ftarget,gpstruct,optimState,0);
                [~,index] = min(z);                
            end
            
            % Something went wrong, random vector
            if isempty(index) || isnan(index); index = randi(size(upoll,1)); end 
            
            % Compute probability that improvement at any location is 
            % less than SufficientImprovement (assumes independence --
            % conservative estimate towards continuing polling)
            gammaz = (gpstruct.ystar - SufficientImprovement - fm)./fs;
            if isfinite(gammaz) & isreal(gammaz)
                fpi = 0.5*erfc(-gammaz/sqrt(2));
                fpi = sort(fpi,'descend');
                pless = prod(1-fpi(1:min(nvars,end)));
            else
                pless = 0;
                UnrelGP = 1;
            end

            % If the probability of improvement at any location is
            % lower than SufficientImprovement, abort iteration IF
            % a good candidate was already found OR if the
            % FASTCONVERGENCE option is on (in this case, a safety flag
            % prevents from doing it twice in a row)                
            if ( UnrelGP && goodpoll ) || ...
               ( ~UnrelGP && ... 
            (goodpoll || ((options.ConsecutiveSkipping || lastskipped < iter-1) && optimState.pollcount >= options.MinFailedPollSteps) ) ...
                    && mean(pless) > 1-options.TolPoI )
%                if (goodpoll || (options.FastConvergence && lastskipped < iter-1) ) ...
%                        && pless > 1-options.TolPoI
                lastskipped = iter;
                break; 
            end

            % Evaluate function and store value
            unew = upoll(index,:);
            [fnew,optimState] = funlogger(funwrapper,unew,optimState,'iter');
            
            % Remove polled vector from set
            upoll(index,:) = [];

            % Save statistics of gp prediction
            optimState.gpstats = ...
                savegpstats(optimState.gpstats,fnew,ymu(:,index),ys(:,index),gpstruct.hypweight);
            
            if options.UncertaintyHandling
                % Add search point to training set
                gpstruct = gpTrainingSet(gpstruct, ...
                    'add', ...
                    unew, ...
                    fnew, ...
                    optimState, ...
                    options, ...
                    0);
                
                % Compute estimated function value at point
                [~,~,fnew,fs2] = gppred(unew,gpstruct);
                fnewsd = sqrt(fs2);
                
                % Computed poll improvement
            else
                fnewsd = 0;                
            end
            
            % Found better objective function value
            if fnew - fpollbest + 0*(fnewsd - fpollbestsd) < 0 
                upollbest = unew;
                fpollbest = fnew;
                fpollhyp = gpstruct.hyp;
                fpollbestsd = fnewsd;
                PollImprovement = fval - fpollbest + 0*(fsd - fpollbestsd);
                if PollImprovement > SufficientImprovement
                    goodpoll = 1;
                end
            end

            optimState.pollcount = optimState.pollcount + 1;

            if optimState.pollcount > nvars*2; break; end
        end % Poll loop
        
        % Evaluate poll
        if (PollImprovement > 0 && options.SloppyImprovement) || ...
                PollImprovement > SufficientImprovement
            polldirection = find(abs(upollbest - ubest) > 1e-12,1); % The sign might be wrong for periodic variables (it's unused anyhow)            
            [ubest,fval,fsd,optimState,gpstruct] = UpdateIncumbent(ubest,fval,fsd,upollbest,fpollbest,fpollbestsd,optimState,gpstruct,options);
            u = ubest;
        end

        if PollImprovement > SufficientImprovement
            % Successful poll
            MeshSizeInteger = min(MeshSizeInteger + 1, options.MaxPollGridNumber);
            SuccessPoll = 1;
            optimState.usuccess = [optimState.usuccess; ubest];
            optimState.fsuccess = [optimState.fsuccess; fval];
        else
            % Failed poll
            if options.AccelerateMesh && iter > options.AccelerateMeshSteps && ...
                fvalhistory(iter-options.AccelerateMeshSteps) - fval < options.TolFun
                MeshSizeInteger = MeshSizeInteger - 2;
            else
                MeshSizeInteger = MeshSizeInteger - 1;
            end            
            
            optimState.SearchSizeInteger = min(optimState.SearchSizeInteger, MeshSizeInteger*options.SearchGridMultiplier - options.SearchGridNumber);
            SuccessPoll = 0;
            
            % Profile plot of iteration
            if strcmpi(options.Plot,'profile') && ~isempty(gpstruct.x)
                % figure(iter);
                hold off;                
                landscapeplot(@(u_) funwrapper(origunits(u_,optimState)), ...
                    u, ...
                    LB, ...
                    UB, ...
                    MeshSize, ...
                    gpstruct, ...
                    [], ...
                    31);
                drawnow;
                if options.Debug
                    display(['GP hyper-parameters at iteration ' num2str(iter) ':']);
                    gpstruct.hyp.mean
                    exp(gpstruct.hyp.cov')
                    if numel(gpstruct.hyp.lik) == 1
                        exp(gpstruct.hyp.lik)
                    else
                        [gpstruct.hyp.lik(1),exp(gpstruct.hyp.lik(2:end))']
                    end
                end
                
                % pause
            end
            
        end

        MeshSize = options.PollMeshMultiplier^MeshSizeInteger;
        optimState.meshsize = MeshSize;
        pollCount = pollCount + 1;

        % Print iteration
        if any(strcmpi(options.Display,{'iter','all'}))
            if SuccessPoll; string = 'Successful poll'; else string = 'Refine grid'; end
            action = [];
            if refitted; if isempty(action); action = 'Train'; else action = [action ', train']; end; end            
            if lastskipped == iter; if isempty(action); action = 'Skip'; else action = [action ', skip']; end; end
            if options.UncertaintyHandling
                fprintf(displayFormat,iter,optimState.funccount,fval,fsd,MeshSize,string,action);                
            else
                fprintf(displayFormat,iter,optimState.funccount,fval,MeshSize,string,action);
            end
        end
        
        % Add polled point to training set
        if ~isempty(unew) && ~options.UncertaintyHandling
            gpstruct = gpTrainingSet(gpstruct, ...
                'add', ...
                unew, ...
                fnew, ...
                optimState, ...
                options, ...
                0);
        end        
        
        %H = fhess(@(xi_) gppred(xi_,gpstruct), gridunits(x,optimState), [], 'step', optimState.searchmeshsize)
        %H2 = fhess(@(x_) funwrapper(x_), x, [], 'step', optimState.searchmeshsize);        
        %[H; H2]
        % [inv(H2); optimState.Binv]
        
    end % Poll step

    % Moved during the poll step
    if goodpoll; gpstruct.post = []; end    
    
    %----------------------------------------------------------------------
    %% Finalize iteration
            
    % Scatter plot of iteration
    if strcmpi(options.Plot,'scatter')
        scatterplot(iter,ubest,fval,action,gpstruct,optimState,options);
    end
        
    if MeshSize < options.TolMesh
        isFinished = 1;
    end
                        
    fvalhistory(iter) = fval;
    fhyp = gpstruct.hyp;
    
    % Check termination conditions
    if optimState.funccount >= options.MaxFunEvals; isFinished = 1; end
    if iter >= options.MaxIter; isFinished = 1; end
    if iter > options.TolStallIters && ...
            fvalhistory(iter-options.TolStallIters) - fval < options.TolFun
        isFinished = 1;
    end
        
    % Store best points at the end of each iteration
    if DoPollStep || isFinished
        optimState.iterList.u(iter,:) = u;
        optimState.iterList.fval(iter,1) = fval;
        optimState.iterList.fsd(iter,1) = fsd;
        optimState.iterList.hyp{iter} = gpstruct.hyp;        
    end
    
    % Re-evaluate all points
    if DoPollStep && options.UncertaintyHandling                        
        if iter > 1
            optimState = reevaluateIterList(optimState,gpstruct,options);

            % Pick best estimate        
            y = optimState.iterList.fval;
            [~,index] = min(y);

            fval = optimState.iterList.fval(index);
            fsd = optimState.iterList.fsd(index);
            u = optimState.iterList.u(index,:);
            fhyp = optimState.iterList.hyp{index};        
        end
    end
    
    if isFinished
        if Restarts > 0 && optimState.funccount < options.MaxFunEvals
            fvalhistory(1:iter-1) = Inf;
            isFinished = 0;
            MeshSizeInteger = 0;
            Restarts = Restarts - 1;
        end
    else
        % Iteration count is increased after the poll step
        if DoPollStep; iter = iter + 1; end        
    end

end

x = origunits(u,optimState);    % Convert back to original space
exitflag = 0;
output.FuncCount = optimState.funccount;
output.IterCount = iter;
output.optimState = optimState;

% Return function evaluation struct (can be reused in future runs)
if nargout > 4
    [~,funValues] = funlogger(funwrapper,u,optimState,'done');
end

% Re-evaluate all best points (skip first iteration)
if options.UncertaintyHandling && iter > 1
    
    optimState = reevaluateIterList(optimState,gpstruct,options);
        
    if 0
        % Order by probability of being minimum
        Nsamples = 1e4;
        y = bsxfun(@plus, optimState.iterList.fval(2:end), ...
            bsxfun(@times,optimState.iterList.fsd(2:end),randn(iter-1,Nsamples)));
        [~,winner] = min(y,[],1);
        winner = winner + 1;

        % Compute exceedance probability (probability of being best value)
        xp = zeros(1,iter);
        for i = 2:iter; xp(i) = sum(winner == i); end
        xp = xp/sum(xp);

        [~,index] = max(xp);        
    else
        % Order by lowest probabilistic upper bound
        y = optimState.iterList.fval + 3*optimState.iterList.fsd;
        [~,index] = min(y);        
    end

    fval = optimState.iterList.fval(index);
    fsd = optimState.iterList.fsd(index);
    u = optimState.iterList.u(index,:);
    x = origunits(u,optimState);
    
    if options.Debug
        optimState.iterList.u
        [optimState.iterList.fval, optimState.iterList.fsd]
    end
    
end

% Return mean and SD of the estimated function value at the optimum
output.fval = fval;
output.fsd = fsd;


end


%--------------------------------------------------------------------------
function gpstats = savegpstats(gpstats,fval,ymu,ys,hypw)
%SAVEGPSTATS Save statistics of GP prediction

nsamples = size(hypw,2);

if isempty(gpstats)
    n = 200;
    gpstats.fval = zeros(1,n);
    gpstats.ymu = zeros(nsamples,n);
    gpstats.ys = zeros(nsamples,n);
    gpstats.hypw = zeros(nsamples,n);
    gpstats.last = 0;
else
    idx = gpstats.last + 1;
    gpstats.fval(idx) = fval;
    gpstats.ymu(1:nsamples,idx) = ymu;
    gpstats.ys(1:nsamples,idx) = ys;
    gpstats.hypw(1:nsamples,idx) = hypw;
    gpstats.last = idx;    
end

end

%--------------------------------------------------------------------------
function [refitgp,UnrelGP,optimState] = IsRefitTime(optimState,options)
%ISREFITTIME Check if it is time to retrain the Gaussian Process

nvars = size(optimState.X,2);

% Check gp prediction statistics (track model reliability)
try
    UnrelGP = gppredcheck(optimState.gpstats,options.NormAlphaLevel);
catch
    UnrelGP = 1;
end

if optimState.funccount < 200
    refitperiod = max(10, nvars*2);
else
    refitperiod = nvars*5;
end

refitgp = optimState.lastfitgp < (optimState.funccount - options.MinRefitTime) && ...
    (optimState.gpstats.last >= refitperiod || UnrelGP) && optimState.funccount > nvars;

if refitgp
    %[plo,t,phi]
    %[gpstruct.hyp.cov(:)',gpstruct.hyp.lik]
    optimState.lastfitgp = optimState.funccount;
    % Reset gp prediction statistics
    optimState.gpstats = savegpstats([],[],[],[],ones(1,max(1,options.gpSamples)));
    UnrelGP = 0;
end

end
              
%--------------------------------------------------------------------------
function [unew,fvalnew,fsdnew,optimState,gpstruct] = UpdateIncumbent(uold,fvalold,fsdold,unew,fvalnew,fsdnew,optimState,gpstruct,options)
%UPDATEINCUMENT Move to a new point

optimState.u = unew;
optimState.fval = fvalnew;
optimState.fsd = fsdnew;

if options.HessianUpdate && ~strcmpi(options.HessianMethod,'cmaes')
    nvars = size(unew,2);
    
    % Compute derivatives at old and new point
    g1 = fgrad(@(x_) gppred(x_,gpstruct), unew,'central','step',optimState.searchmeshsize);
    g0 = fgrad(@(x_) gppred(x_,gpstruct), uold,'central','step',optimState.searchmeshsize);
    g1 = g1./optimState.scale; g0 = g0./optimState.scale;
    optimState.grad = g1(:); % Gradient at new incumbent
    y = g1 - g0;
    s = (unew - uold)./optimState.scale;
    % sqrt(sum(s.*s))
    
    if any(optimState.periodicvars)
        % Compute vector difference as the closest distance along periodic dimensions    
        % dual = optimState.UB - optimState.LB - abs(s);   % Opposite distance
        dual = (optimState.UB - optimState.LB)./optimState.scale - abs(s);   % Opposite distance
        for dim = find(optimState.periodicvars)
            if dual(dim) < abs(s(dim))      % Opposite distance is shorter
                s(dim) = -sign(s(dim))*dual(dim);
            end
        end
    end
    
    Binv = optimState.Binv;
    switch options.HessianMethod
        case {'bfsg','bfgs','bgfs','bgsf','bsfg','bsgf'} % I never get it right
            if y*s' > 0
                Binv = Binv + (s*y' + y*Binv*y')*(s'*s)/(s*y')^2 - (Binv*y'*s + s'*y*Binv)/(s*y');
            else
                optimState.violations = optimState.violations + 1;
                % if optimState.violations < 1; return; end
                display(['negative ' num2str(optimState.violations)]);
                Binv = diag(gpstruct.lenscale.^2);
                optimState.violations = 0;
            end
            % Binv = 0.95*Binv + 0.05*eye(nvars) + (s*y' + y*Binv*y')*(s'*s)/(s*y')^2 - (Binv*y'*s + s'*y*Binv)/(s*y');
        case 'naive'
            c1 = 2/nvars^2; c2 = 0; c0 = c1+c2; c0 = 0;
            % gn = g1/sqrt(g1*g1');
            % Binv = Binv*(1-c0) + c1*(s'*s)/optimState.meshsize^2;
            Binv = Binv*(1-c0) + c1*(s'*s)/(s*s');
            if c2 > 0
                % Binv = Binv + c2*(gn'*gn)*((s*s')/optimState.meshsize)^2;            
                Binv = Binv + c2*diag(gpstruct.pollscale.^2); % /sqrt(mean(gpstruct.pollscale.^2));            
            end
            % s/optimState.meshsize
            % gn*optimState.meshsize
        case 'hessian'
            B = fhess(@(xi_) gppred(xi_,gpstruct), unew, [], 'step', optimState.searchmeshsize);
            c1 = 0;
            Binv = (1-c1)*Binv + c1*inv(B);
        case 'neighborhood'
            
            if 1
                optionstemp = options;
                optionstemp.Ndata = 2^10;
                optionstemp.MinNdata = 4 + floor(3*log(nvars));
                optionstemp.BufferNdata = Inf;
                optionstemp.gpRadius = 2;
                gptemp = gpTrainingSet(gpstruct, ...
                    'neighborhood', ...
                    uold, ...
                    [], ...
                    optimState, ...
                    optionstemp, ...
                    0);
                muratio = 0.25;
            else
                gptemp = gpstruct;
                muratio = 0.5;
            end
            X = gptemp.x;
            Y = gptemp.y;

            % Small jitter added to each direction
            jit = optimState.searchmeshsize;

            % Compute vector weights
            mu = floor(muratio*size(X,1));
            mu = max([mu, (4 + floor(3*log(nvars)))/2, floor(size(X,1)/2)]);
            weights = zeros(1,1,floor(mu));
            weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
            weights = weights./sum(weights);

            % Compute top vectors
            [~,index] = sort(Y,'ascend');

            % Compute weighted covariance matrix wrt X0
            u = uold;
            Ubest = X(index(1:floor(mu)),:);
            C = ucov(Ubest,u,weights,optimState)./optimState.meshsize^2;
                        
            mueff = 1/sum(weights.^2);
            amu = 2; c1 = 0;
            % c1 = 2/((nvars+1.3)^2 + mueff); % Doesn't seem to improve --
            % might try to implement the whole path thing
            cmu = min(1-c1, amu*(mueff-2+1/mueff)/((nvars+2)^2 + amu*mueff/2));
                        
            Binv = (1-cmu-c1)*Binv + c1*(s'*s)/optimState.meshsize^2 + cmu*C;
            
        otherwise
            error('Unknown Hessian update method.');
    end
    
    optimState.Binv = Binv;
    
    if mod(optimState.iter,2) == 0 && options.HessianAlternate
        optimState.C = diag(gpstruct.pollscale./sqrt(sum(gpstruct.pollscale.^2)));
    else        
        try
            [V,D] = eig(Binv);
            lambda = real(diag(D));
            % lambda(:)'
        catch
            lambda = [];
        end
        if isempty(lambda) || (min(lambda) < 0 && abs(min(lambda)/max(lambda)) > 1e-14)
            optimState.Binv = eye(nvars);
            optimState.B = eye(nvars);
            optimState.C = eye(nvars);
            display('reset')
        else
            lambda = max(lambda, max(lambda)*1e-14);
            optimState.B = real(V)*diag(1./lambda)*real(V)';
            lambdared = sqrt(lambda/sum(lambda));   
            lambdared = min(max(lambdared, optimState.searchmeshsize), max((optimState.UB-optimState.LB)./optimState.scale));
            % optimState.B = real(V)*diag(1./lambdared.^2)*real(V)';
            lambdared = lambdared/sqrt(sum(lambdared.^2));
            optimState.C = real(V)*diag(lambdared);
            optimState.Cres = optimState.C/sqrt(sum(1./lambdared));
                        
            if any(~isreal(optimState.C(:)))
                optimState.Binv = eye(nvars);
                optimState.C = eye(nvars);
                display('unreal')
            end
            
            % log10(lambda)
        end
    end    
end

end
            
%UPDATEINCUMBENTPREDICTION Update gp mean and variance at incumbent
function optimState = UpdateIncumbentPrediction(ubest,hyp,optimState,gpstruct,options)

if options.UncertaintyHandling || options.UncertainIncumbent
    gptemp = gpstruct;
    gptemp.hyp = hyp;
    [~,~,ftargetmu,ftargets2] = gppred(ubest,gptemp);
    ftargetmu = real(ftargetmu);
    ftargets = sqrt(max(ftargets2,0));
    if ~isfinite(ftargetmu) || ~isreal(ftargets) || ~isfinite(ftargets)
        ftargetmu = optimState.fval;
        ftargets = optimState.fsd;
    end
    
    if ~options.AcqHedge
        if options.AlternativeIncumbent
            % [ftargetmu - optimState.sdlevel*ftargets - optimState.fval]            
            % ftarget = max(ftargetmu - optimState.sdlevel*ftargets, optimState.fval) - options.TolFun;                    
            ftarget = ftargetmu - sqrt(size(ubest,2))/sqrt(optimState.funccount)*ftargets;
        else
            ftarget = ftargetmu - optimState.sdlevel*sqrt(ftargets2 + options.TolFun^2);
        end
    end
    % ftarget = ftarget - log(optimState.funccount)*sqrt(fs2);
else
    ftarget = optimState.fval - options.TolFun;
    ftargetmu = optimState.fval;
    ftargets = 0;
end

optimState.ftargetmu = ftargetmu;
optimState.ftargets = ftargets;        
optimState.ftarget = ftarget;

end


%--------------------------------------------------------------------------
function optimState = UpdateSearch(optimState,searchstatus,searchdist,options)
%UPDATESEARCH Update statistics of searches

% Initialize statistics
if ~isfield(optimState,'searchstats') || isempty(optimState.searchstats)
    optimState.searchstats.logsearchfactor = [];
    optimState.searchstats.success = [];
    optimState.searchstats.udist = [];
end

optimState.searchstats.logsearchfactor(end+1) = log(optimState.searchfactor);
optimState.searchstats.udist(end+1) = searchdist;

switch searchstatus
    case 'success'
        optimState.searchstats.success(end+1) = 1;
        optimState.searchfactor = optimState.searchfactor*options.SearchScaleSuccess;
        if options.AdaptiveIncumbentShift; optimState.sdlevel = optimState.sdlevel*2; end
    case 'incremental'
        optimState.searchstats.success(end+1) = 0.5;
        optimState.searchfactor = optimState.searchfactor*options.SearchScaleIncremental;
        if options.AdaptiveIncumbentShift; optimState.sdlevel = optimState.sdlevel*(2^2); end
    case 'failure'
        optimState.searchstats.success(end+1) = 0;
        optimState.searchfactor = max(options.SearchFactorMin,optimState.searchfactor*options.SearchScaleFailure);
        if options.AdaptiveIncumbentShift; optimState.sdlevel = max(options.IncumbentSigmaMultiplier,optimState.sdlevel/2); end
end

% optimState.sdlevel

% Reset search factor at the end of each search
if optimState.searchcount == options.SearchNtry
    index = max(1,numel(optimState.searchstats.success)-100):numel(optimState.searchstats.success);
    success = optimState.searchstats.success(index);
    ud = optimState.searchstats.udist(index);
    if sum(success == 1) > 3 && ...
            sum(success == 0.5) > 1 && ...
            sum(success == 0) > 3 && 0
        
        msucc = median(ud(success == 1));
        mincr = prctile(ud(success == 0.5),90);
        mfail = median(ud(success == 0));
        
        [msucc mincr mfail]
        
        %optimState.searchfactor = ...
        %    exp(0.25*log(msucc) + 0.25*log(mincr) + 0.5*log(mfail));
        %if optimState.searchfactor < 2*mincr
            optimState.searchfactor = ...
                4*max(exp(0.5*log(msucc) + 0.5*log(mfail)), 4*mincr);
        %end
        
        optimState.searchfactor
    else
        optimState.searchfactor = 1;        
    end
end

% optimState.searchfactor = exp(randn());

end

%--------------------------------------------------------------------------
function optimState = reevaluateIterList(optimState,gpstruct,options)
%REEVALUATEITERLIST Re-evaluate function values at stored locations

iter = optimState.iter;

for index = 1:iter
    gpstruct.hyp = optimState.iterList.hyp{index};

    ui = optimState.iterList.u(index,:);

    gpstruct = gpTrainingSet(gpstruct, ...
        options.gpMethod, ...
        ui, ...
        [], ...
        optimState, ...
        options, ...
        0);

    % Compute estimated function value at point
    [~,~,fval,fs2] = gppred(ui,gpstruct);
    fsd = sqrt(fs2);

    optimState.iterList.fval(index) = fval;
    optimState.iterList.fsd(index) = fsd;
end

end

%--------------------------------------------------------------------------
function optimState = updateSearchBounds(optimState)
%UPDATESEARCHBOUNDS Update bounds for search on the mesh.

optimState.LBsearch = force2grid(optimState.LB,optimState);
optimState.LBsearch(optimState.LBsearch < optimState.LB) = ...
    optimState.LBsearch(optimState.LBsearch < optimState.LB) + optimState.searchmeshsize;
optimState.UBsearch = force2grid(optimState.UB,optimState);
optimState.UBsearch(optimState.UBsearch > optimState.UB) = ...
    optimState.UBsearch(optimState.UBsearch > optimState.UB) - optimState.searchmeshsize;
end

%--------------------------------------------------------------------------
function [hedge, acqIndex, ymu, ys, fmu, fs] = acqPortfolio(method, hedge, u, f, fs, gpstruct, optimState, options, SufficientImprovement, fvalold, MeshSize)

if nargin < 9; fvalold = []; end
if nargin < 10; MeshSize = []; end

switch lower(method(1:3))
    case 'acq'
        
        [acqIndex,~,ymu,ys,fmu,fs] = ...
            acqHedge(u,f,fs,gpstruct,optimState,options,SufficientImprovement);

        % Initialize hedge struct
        if isempty(hedge)
            [hedge.g,hedge.str] = acqHedge();
            hedge.n = numel(hedge.g);
            hedge.count = 0;
            hedge.gamma = options.HedgeGamma;
            hedge.beta = options.HedgeBeta;
            hedge.decay = options.HedgeDecay;
            hedge.string = {'mpi','mei'};
        end

        hedge.count = hedge.count + 1;
        % gammaHedge = min(1/nHedge, sqrt(log(nHedge)/(nHedge*hedge.count)));

        hedge.p = exp(hedge.beta*(hedge.g - max(hedge.g)))./sum(exp(hedge.beta*(hedge.g - max(hedge.g))));
        hedge.p = hedge.p*(1-hedge.n*hedge.gamma) + hedge.gamma;
        % hedge.p
        hedge.chosen = find(rand() < cumsum(hedge.p),1);
        if hedge.gamma == 0
            hedge.phat = ones(size(hedge.g));
        else
            hedge.phat = Inf(size(hedge.p));
            hedge.phat(hedge.chosen) = hedge.p(hedge.chosen);
        end
        
    case 'upd'

        for iHedge = 1:hedge.n
            uHedge = u(min(iHedge,end),:);

            %gpstructnew = gpTrainingSet(gpstructnew, ...
            %    options.gpMethod, ...
            %    uHedge, ...
            %    [], ...
            %    optimState, ...
            %    options, ...
            %    0);

            if iHedge == hedge.chosen
                fHedge = f; fsHedge = fs;
            elseif hedge.gamma == 0
                % Compute estimated function value at point
                [~,~,fHedge,fs2Hedge] = gppred(uHedge,gpstruct);
                fsHedge = sqrt(fs2Hedge);                
            else
                fHedge = 0; fsHedge = 1;
            end
            
            if fsHedge == 0
                er = max(0, fvalold - fHedge);
            elseif isfinite(fHedge) && isfinite(fsHedge) && isreal(fsHedge) && fsHedge > 0
                % Probability of improvement
                gammaz = (fvalold - fHedge)./fsHedge;
                fpi = 0.5*erfc(-gammaz/sqrt(2));            

                % Expected reward
                er = fsHedge.*(gammaz.*fpi + exp(-0.5*(gammaz.^2))/sqrt(2*pi));
            else
                er = 0;
            end

            hedge.g(iHedge) = hedge.decay*hedge.g(iHedge) + er/hedge.phat(iHedge)/MeshSize;
        end
         
end
end
