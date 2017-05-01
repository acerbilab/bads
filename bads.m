function [x,fval,exitflag,output,optimState,gpstruct] = bads(fun,x0,LB,UB,PLB,PUB,nonbcon,options,varargin)
%BADS Constrained optimization using Bayesian Adaptive Direct Search
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
%      - if X0 is empty, LB and UB need to be specified as vectors.
%
%   X = BADS(FUN,X0,LB,UB,PLB,PUB) specifies a set of plausible lower and
%   upper bounds such that LB <= PLB <= X0 <= PUB <= UB. Both PLB and PUB
%   need to be finite. PLB and PUB are used to design the initial mesh of 
%   the direct search, and represent a plausible range for the 
%   optimization variables. As a rule of thumb, set PLB and PUB such that 
%   there is > 90% probability that the minimum is found within the box 
%   (where in doubt, just set PLB=LB and PUB=UB).
%
%   X = BADS(FUN,X0,LB,UB,PLB,PUB,NONBCON) subjects the minimization to the 
%   non-bound constraints defined in NONBCON. The function NONBCON accepts 
%   a N-by-D matrix XI where N is any number of points to evaluate and D is
%   the number of dimensions, and returns the vector C, representing the
%   degree of violation of non-bound inequalities for each point in XI. 
%   BADS minimizes FUN such that C(X)<=0.
%
%   X = BADS(FUN,X0,LB,UB,PLB,PUB,NONBCON,options) minimizes with the default 
%   optimization parameters replaced by values in the structure OPTIONS.
%   BADS('defaults') returns the default OPTIONS struct.
%
%   [X,FVAL] = BADS(...) returns FVAL, the value of the objective function 
%   FUN at the solution X. If the target function is stochastic, FVAL is
%   the expected mean of the function value at X.
%
%   [X,FVAL,EXITFLAG] = BADS(...) returns EXITFLAG which describes the exit 
%   condition of BADS. Possible values of EXITFLAG and the corresponding 
%   exit conditions are:
%
%     0  Maximum number of function evaluations or iterations reached.
%     1  Magnitude of mesh size is less than specified tolerance.
%     2  Change in estimated function value less than the specified tolerance 
%        for OPTIONS.TolStallIters iterations.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = BADS(...) returns a structure OUTPUT with the
%   following information:
%          function: <Objective function name>
%       problemtype: <Type of problem> (unconstrained or bound constrained)
%        targettype: <Type of target function> (deterministic or stochastic)
%        iterations: <Total iterations>
%         funccount: <Total function evaluations>
%          meshsize: <Mesh size at X>
%          overhead: <Fractional overhead (total runtime / total fcn time - 1)>
%          rngstate: <Status of random number generator>
%         algorithm: <Bayesian adaptive direct search>
%           message: <BADS termination message>
%              fval: <Expected mean of function value at X>
%               fsd: <Estimated standard deviation of function value at X>
%
%   [X,FVAL,EXITFLAG,OUTPUT,OPTIMSTATE] = BADS(...) returns a detailed
%   optimization structure OPTIMSTATE.
%
%   [X,FVAL,EXITFLAG,OUTPUT,OPTIMSTATE,GPSTRUCT] = BADS(...) returns the
%   Gaussian Process structure GPSTRUCT.
%
%   OPTIONS = BADS('defaults') returns a basic default OPTIONS structure.
%
%   Optimization of noisy functions:
%   To run BADS on a noisy (stochastic) objective function, set
%       OPTIONS.UncertaintyHandling = 1
%       OPTIONS.NoiseSize = SIGMA
%   where SIGMA is an estimate of the SD of the noise in your problem in
%   a good region of the parameter space. (If not specified, default 
%   SIGMA = 1).
%   Set OPTIONS.UncertaintyHandling = 0 for a deterministic function.
%   If OPTIONS.UncertaintyHandling is not specified, BADS will determine at
%   runtime if the objective function is noisy.
%
%   THIS VERSION IS FOR INTERNAL USE ONLY -- PLEASE DO NOT REDISTRIBUTE.
%   Contact <luigi.acerbi@nyu.edu> for more information.

%--------------------------------------------------------------------------
% BADS: Bayesian Adaptive Direct Search for nonlinear function minimization.
% To be used under the terms of the GNU General Public License 
% (http://www.gnu.org/copyleft/gpl.html).
%
%   Author (copyright): Luigi Acerbi, 2017
%   e-mail: luigi.acerbi@{gmail.com,nyu.edu}
%   URL: http://luigiacerbi.com
%   Release date: Apr 29, 2017
%   Version: 0.95
%   References: Check https://github.com/lacerbi/bads
%--------------------------------------------------------------------------

% TO DO:
% - check options.FunValues provided as argument
% - play around with OPTIONS.ImprovementQuantile
% - added retrain at first search step (is it necessary?)
% - removed added point at the end of poll stage (was it necessary?)
% - fix restarts e multibayes
% - compute func running time and do more stuff if func is slow
% - understand if the warping is working correctly, test moar
% - aggiungi warning se la mesh cerca di espandere oltre MaxPollGridNumber
%   (sintomo di misspecification dei bound)

%% Start timer

t0 = tic;

%% Basic default options

defopts.Display                 = 'iter         % Level of display ("iter", "notify", "final", or "off")';
defopts.MaxIter                 = '200*nvars    % Max number of iterations';
defopts.MaxFunEvals             = '500*nvars    % Max number of objective fcn evaluations';
defopts.FunValues               = '[]           % Struct with pregress fcn evaluations (X and Y fields)';
defopts.PeriodicVars            = '[]           % Array with indices of periodic variables';
defopts.NonlinearScaling        = 'on           % Automatic nonlinear rescaling of variables';
defopts.CompletePoll            = 'off          % Complete polling around the current iterate';
defopts.AccelerateMesh          = 'on           % Accelerate mesh contraction';
defopts.UncertaintyHandling     = '[]           % Explicit noise handling (if empty, determine at runtime)';
defopts.NoiseObj                = 'off          % Objective fcn returns noise estimate as 2nd argument (unsupported)';
defopts.NoiseSize               = '[]           % Base observation noise magnitude';
defopts.OptimToolbox            = '[]           % Use Optimization Toolbox (if empty, determine at runtime)';

%% If called with no arguments or with 'defaults', return default options
if nargin < 1 || strcmpi(fun,'defaults')
    if nargin < 1
        fprintf('Basic default options returned (type "help bads" for help).\n');
    end
    x = defopts;
    return;
end

%% Advanced options (do not modify unless you *know* what you are doing)

defopts.Plot                    = 'off                  % Show optimization plots ("profile", "scatter", or "off")';
defopts.Debug                   = 'off                  % Debug mode, plot additional info';

% Termination conditions
defopts.TolMesh                 = '1e-6                 % Tolerance on mesh size';
defopts.TolFun                  = '1e-3                 % Min significant change of objective fcn';
defopts.TolStallIters           = '4 + floor(nvars/2)   % Max iterations with no significant change (doubled under uncertainty)';
defopts.TolNoise                = 'sqrt(eps)*options.TolFun  % Min variability for a fcn to be considered noisy';

defopts.Ninit                   = 'nvars                % Number of initial objective fcn evaluations';
defopts.InitFcn                 = '@initSobol           % Initialization function';
% defoptions.InitFcn            = '@initLHS';
defopts.PollMethod              = '@pollMADS2N          % Poll function';
defopts.Nbasis                  = '200*nvars';

defopts.Restarts                = '0                    % Number of restart attempts';
defopts.CacheSize               = '1e4                  % Size of cache for storing function evaluations';

defopts.TolImprovement          = '1                    % Minimum significant improvement at unit mesh size';
defopts.ForcingExponent         = '3/2                  % Exponent of forcing function';
defopts.PollMeshMultiplier      = '2                    % Mesh multiplicative factor between iterations';
defopts.IncumbentSigmaMultiplier = '0.1                 % Multiplier to incumbent uncertainty for acquisition functions';
defopts.ImprovementQuantile     = '0.5                  % Quantile when computing improvement (<0.5 for conservative improvement)';
defopts.FinalQuantile           = '1e-3                 % Top quantile when choosing final iteration';
defopts.ForcePollMesh           = 'no                   % Force poll vectors to be on mesh';

defopts.AlternativeIncumbent    = 'off                  % Use alternative incumbent offset';
defopts.AdaptiveIncumbentShift  = 'off                  % Adaptive multiplier to incumbent uncertainty';
defopts.gpRescalePoll           = '1                    % GP-based geometric scaling factor of poll vectors';
defopts.FitnessShaping          = 'off                  % Nonlinear rescaling of objective fcn';
defopts.WarpFunc                = '0                    % GP warping function type';

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
defopts.MinNdata                = '50                   % Minimum number of training data (doubled under uncertainty)';
defopts.BufferNdata             = '100                  % Max number of training data removed if too far from current point';
defopts.gpSamples               = '0                    % Hyperparameters samples (0 = optimize)';
defopts.MinRefitTime            = '2*nvars              % Minimum fcn evals before refitting the GP';
defopts.PollTraining            = 'yes                  % Train GP also during poll stage';
defopts.DoubleRefit             = 'off                  % Always try a second GP fit';
defopts.gpMeanPercentile        = '90                   % Percentile of empirical GP mean';
defopts.gpMeanRangeFun          = '@(ym,y) (ym - prctile(y,50))/5*2   % Empirical range of hyperprior over the mean';
defopts.gpdefFcn                = '{@gpdefBads,''rq'',[1,1]}  % GP definition fcn';
defopts.gpMethod                = 'nearest              % GP training set selection method';
% defopts.gpMethod                = 'grid                 % GP training set selection method';
defopts.gpCluster               = 'no                   % Cluster additional points during training';
defopts.RotateGP                = 'no                   % Rotate GP basis';
defopts.gpRadius                = '3                    % Radius of training set';
defopts.UseEffectiveRadius      = 'yes                   %';
defopts.gpCovPrior              = 'iso                  % GP hyper-prior over covariance';
defopts.gpFixedMean             = 'no';
defopts.FitLik                  = 'yes                  % Fit the likelihood term';
defopts.PollAcqFcn              = '{@acqLCB,[]}         % Acquisition fcn for poll stage';
defopts.SearchAcqFcn            = '{@acqLCB,[]}         % Acquisition fcn for search stage';
defopts.AcqHedge                = 'off                  % Hedge acquisition function';
defopts.CholAttempts            = '0                    % Attempts at performing the Cholesky decomposition';
defopts.NoiseNudge              = '[1 0]                % Increase nudge to noise in case of Cholesky failure';
defopts.RemovePointsAfterTries  = '1                    % Start removing training points after this number of failures';
defopts.gpSVGDiters             = '200                  % SVGD iterations for GP training';
defopts.gpWarnings              = 'off          % Issue warning if GP hyperparameters fit fails';

defopts.UncertainIncumbent      = 'yes                  % Treat incumbent as if uncertain regardless of uncertainty handling';
defopts.MeshNoiseMultiplier     = '0.5                  % Contribution to log noise magnitude from log mesh size (0 for noisy functions)';
defopts.TolPoI                  = '1e-6/nvars           % Threshold probability of improvement (PoI); set to 0 to always complete polling';
defopts.SkipPoll                = 'yes                  % Skip polling if PoI below threshold, even with no success';
defopts.ConsecutiveSkipping     = 'yes                  % Allow consecutive incomplete polls';
defopts.SkipPollAfterSearch     = 'yes                  % Skip polling after successful search';
defopts.MinFailedPollSteps      = 'Inf                  % Number of failed fcn evaluations before skipping is allowed';
defopts.NormAlphaLevel          = '1e-6                 % Alpha level for normality test of gp predictions';
defopts.AccelerateMeshSteps     = '3                    % Accelerate mesh after this number of stalled iterations';
defopts.SloppyImprovement       = 'yes                  % Move incumbent even after insufficient improvement';
defopts.HessianUpdate           = 'no                   % Update Hessian as you go';
defopts.HessianMethod           = 'bfgs                 % Hessian update method';
defopts.HessianAlternate        = 'no                   % Alternate Hessian iterations';

defopts.HedgeGamma              = '0.125';
defopts.HedgeBeta               = '1e-3/options.TolFun';
defopts.HedgeDecay              = '0.1^(1/(2*nvars))';

defopts.TrueMinX                = '[]                   % Location of the global minimum (for visualization only)';


%% If called with 'all', return all default options
if strcmpi(fun,'all')
    x = defopts;
    return;
end

%% Check that all BADS subfolders are on the MATLAB path
add2path();

%% Input arguments

if nargin < 3 || isempty(LB); LB = -Inf; end
if nargin < 4 || isempty(UB); UB = Inf; end
if nargin < 5; PLB = []; end
if nargin < 6; PUB = []; end
if nargin < 7; nonbcon = []; end
if nargin < 8; options = []; end

% Retro-compatibility with older interface without NONBCON
if ~isempty(nonbcon) && isstruct(nonbcon)
    if ~isempty(options)
        nvarargin = numel(varargin);
        for i = nvarargin+1:-1:2; varargin{i} = varargin{i-1}; end
        varargin{1} = options;
    end
    options = nonbcon;
    nonbcon = [];
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
end

nvars = numel(x0);
optimState = [];

% Setup algorithm options
options = setupoptions(nvars,defopts,options);    

% Setup and transform variables
[u0,LB,UB,PLB,PUB,MeshSizeInteger,optimState] = ...
    setupvars(x0,LB,UB,PLB,PUB,optimState,options,prnt);
    
optimState = updateSearchBounds(optimState);

% Store objective function
if ischar(fun); fun = str2func(fun); end
optimState.fun = fun;
if isempty(varargin)
    funwrapper = fun;   % No additional function arguments passed
else
    funwrapper = @(u_) fun(u_,varargin{:});
end

% Store constraints function
if ischar(nonbcon); nonbcon = str2func(nonbcon); end
optimState.nonbcon = nonbcon;

% Initialize function logger
[~,optimState] = funlogger([],u0,optimState,'init',options.CacheSize,options.NoiseObj);

%% Initial function evaluations

iter = 0;
optimState.iter = iter;

% Evaluate starting point and initial mesh, determine if function is noisy
[u,fval,isFinished_flag,optimState,displayFormat] = ...
    evalinitmesh(u0,funwrapper,optimState,options,prnt);
if ~isfinite(fval); error('Cannot find valid starting point.'); end
exitflag = 0;
msg = 'Optimization terminated: reached maximum number of function evaluations after initialization.';
    
% Change options for uncertainty handling
if optimState.UncertaintyHandling
    options.TolStallIters = 2*options.TolStallIters;
    options.Ndata = max(200,options.Ndata);
    options.MinNdata = 2*options.MinNdata;
    %options.gpMeanPercentile = 50;
    options.MinFailedPollSteps = Inf;
    options.MeshNoiseMultiplier = 0;
    if isempty(options.NoiseSize); options.NoiseSize = 1; end
else
    if isempty(options.NoiseSize); options.NoiseSize = sqrt(options.TolFun); end
end

if optimState.UncertaintyHandling  % Current uncertainty in estimate
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
gpstruct = feval(options.gpdefFcn{:},nvars,gplik,optimState,options,[]);
gpstruct.fun = funwrapper;
fhyp = gpstruct.hyp;

% Initialize struct with GP prediction statistics
Nsamples = max(1,options.gpSamples);
optimState.gpstats = savegpstats([],[],[],[],ones(1,Nsamples)/Nsamples);

lastskipped = 0;                % Last skipped iteration

SearchSuccesses = 0;
SearchSpree = 0;
Restarts = options.Restarts;

%% Optimization loop
iter = 1;

while ~isFinished_flag
    optimState.iter = iter;
    refitted_flag = false;  % GP refitted this iteration
    gpexitflag = Inf;       % Exit flag from GP training
    action = [];            % Action performed this iteration (for printing purposes)
                
    % Compute mesh size and search mesh size
    MeshSize = options.PollMeshMultiplier^MeshSizeInteger;
    if options.SearchSizeLocked
        optimState.SearchSizeInteger = min(0,MeshSizeInteger*options.SearchGridMultiplier - options.SearchGridNumber);
    end
    optimState.meshsize = MeshSize;
    optimState.searchmeshsize = options.PollMeshMultiplier.^optimState.SearchSizeInteger;
    
    % Update bounds to grid for search mesh
    optimState = updateSearchBounds(optimState);
    
    % Minimum improvement for a poll/search to be considered successful
    SufficientImprovement = options.TolImprovement*(MeshSize^options.ForcingExponent);
    if options.SloppyImprovement
        SufficientImprovement = max(SufficientImprovement, options.TolFun);
    end

    optimState.SearchSufficientImprovement = SufficientImprovement;
    % Multiple successful searches raise the bar for improvement
    % optimState.SearchSufficientImprovement = SufficientImprovement*(2^SearchSuccesses);
        
    %----------------------------------------------------------------------
    %% Search stage
            
    % Perform search if there are still available attempts, and if there 
    % are more than NVARS stored points
    DoSearchStep_flag = optimState.searchcount < options.SearchNtry ...
        && size(gpstruct.y,1) > nvars;
    
    if DoSearchStep_flag
        
        % Check whether it is time to refit the GP
        [refitgp_flag,~,optimState] = IsRefitTime(optimState,options);
        if refitgp_flag || optimState.searchcount == 0; gpstruct.post = []; end
        
        if isempty(gpstruct.post)
            % Local GP approximation on current point
            [gpstruct,gptempflag] = gpTrainingSet(gpstruct, ...
                options.gpMethod, ...
                u, ...
                [], ...    %             [upoll; gridunits(x,optimState)], ...
                optimState, ...
                options, ...
                refitgp_flag);
                if refitgp_flag; refitted_flag = true; end
                gpexitflag = min(gptempflag,gpexitflag);
        end
        
        % Update optimization target (based on GP prediction at incumbent)
        optimState = UpdateTarget(ubest,fhyp,optimState,gpstruct,options);        
        
        % Generate search set (normalized coordinates)        
        optimState.searchcount = optimState.searchcount + 1;        
        [usearchset,optimState] = feval(options.SearchMethod{:}, ...
            u, ...
            gpstruct, ...
            LB, ...
            UB, ...
            optimState, ...
            options);    
        
        % Enforce periodicity
        usearchset = periodCheck(usearchset,LB,UB,optimState);
        
        % Force candidate points on search grid
        usearchset = force2grid(usearchset, optimState);
                
        % Remove already evaluated or unfeasible points from search set
        usearchset = uCheck(usearchset,optimState.TolMesh,optimState,1);
        
        if ~isempty(usearchset) % Non-empty search set

            ns = size(usearchset, 1);
            ymu = zeros(numel(gpstruct.hyp),ns);
            ys = zeros(numel(gpstruct.hyp),ns);
            
            % Evaluate acquisition function on search set
            try
                %----------------------------------------------------------
                if options.AcqHedge
                    [optimState.hedge,acqIndex,ymu,ys] = ...
                        acqPortfolio('acq',optimState.hedge,usearchset,optimState.ftarget,0,gpstruct,optimState,options,SufficientImprovement);
                    index = acqIndex(optimState.hedge.chosen);
                    z = 1;
                %----------------------------------------------------------
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
                        
            acqu = [];
            
            %--------------------------------------------------------------
            % Local optimization of the acquisition function 
            % (generally it does not improve results)
            if options.SearchOptimize
                acqoptoptions = optimset('Display','off','GradObj','off','DerivativeCheck','off',...
                    'TolX',optimState.TolMesh,'TolFun',options.TolFun);

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
                    acqu = uCheck(acqu,optimState.TolMesh,optimState,1);
                catch
                    acqu = [];
                end
            end
            %--------------------------------------------------------------
            
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
            
            if optimState.UncertaintyHandling
                gpstructnew = gpTrainingSet(gpstruct, ...
                    options.gpMethod, ...
                    usearch, ...
                    [], ...
                    optimState, ...
                    options, ...
                    0);
                
                % Compute estimated function value at point
                [~,~,fsearch,fs2] = gppred(usearch,gpstructnew);                                
                if numel(gpstructnew.hyp) > 1
                    fsearch = weightedsum(gpstructnew.hypweight,fsearch,1);
                    fs2 = weightedsum(gpstructnew.hypweight,fs2,1);
                end
                fsearchsd = sqrt(fs2);
            else
                fsearchsd = 0;
            end
            
            % Compute distance of search point from current point
            searchdist = sqrt(udist(ubest,usearch,gpstruct.lenscale,optimState));            
            
        else    % Empty search set
            fsearch = fval;
            fsearchsd = 0;
            searchdist = 0;
        end
        
        %------------------------------------------------------------------
        % CMA-ES like estimation of local covariance structure (unused)
        if options.HessianUpdate && strcmpi(options.HessianMethod,'cmaes')
            optimState = covmatadapt(u,LB,UB,gpstruct,optimState,options);
        end
        %------------------------------------------------------------------
                
        % Evaluate search
        
        SearchImprovement = EvalImprovement(fval,fsearch,fsd,fsearchsd,options.ImprovementQuantile);
        fvalold = fval;
        
        % Declare if search was success or failure
        if (SearchImprovement > 0 && options.SloppyImprovement) ...
                || SearchImprovement > optimState.SearchSufficientImprovement
            % Search did not fail
            %--------------------------------------------------------------
            if options.AcqHedge
                method = optimState.hedge.str{optimState.hedge.chosen};
            else
            %--------------------------------------------------------------                
                method = feval(options.SearchMethod{:},[],[],[],[],optimState);
            end
            if SearchImprovement > optimState.SearchSufficientImprovement
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
            
            if optimState.UncertaintyHandling; gpstruct = gpstructnew; end
            gpstruct.post = []; % Reset posterior
        else
            % Search failed
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
                        
        %------------------------------------------------------------------
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
        %------------------------------------------------------------------
                
        % Update search statistics and search scale factor
       optimState = UpdateSearch(optimState,searchstatus,searchdist,options);
       
       % Print search results       
       if prnt > 2 && ~isempty(searchstring)
           if optimState.UncertaintyHandling
               fprintf(displayFormat,iter,optimState.funccount,fval,fsd,MeshSize,searchstring,'');                
           else
               fprintf(displayFormat,iter,optimState.funccount,fval,MeshSize,searchstring,'');
           end
       end
       
    end % Search stage

    % Decide whether to perform the poll stage
    switch optimState.searchcount
        case {0, options.SearchNtry}    % Skipped or just finished search             
            optimState.searchcount = 0;
            if SearchSuccesses > 0 && options.SkipPollAfterSearch
                DoPollStep_flag = false;
                SearchSpree = SearchSpree + 1;
                if options.SearchMeshExpand > 0 && ...
                        mod(SearchSpree,options.SearchMeshExpand) == 0
                    MeshSizeInteger = min(MeshSizeInteger + options.SearchMeshIncrement, options.MaxPollGridNumber);
                end
            else
                DoPollStep_flag = true;
                SearchSpree = 0;
            end
            SearchSuccesses = 0;
            % optimState.searchfactor = 1;
        otherwise                       % In-between searches, no poll
            DoPollStep_flag = false;
    end
            
    %----------------------------------------------------------------------
    %% Poll stage  

    u = ubest;
    
    if DoPollStep_flag
        
        PollBestImprovement = 0;        % Best improvement so far
        upollbest = u;                  % Best poll point
        fpollbest = fval;               % Objective value at best point
        fpollhyp = fhyp;                % gp hyper-parameters at best point
        fpollbestsd = fsd;              % Uncertainty of objective func
        optimState.pollcount = 0;       % Poll iterations
        goodpoll_flag = false;          % Found a good poll
        B = [];                         % Poll basis
        upoll = [];                     % Poll vectors
        unew = [];

        % Poll loop
        while (~isempty(upoll) || isempty(B)) ...
                && optimState.funccount < options.MaxFunEvals ...
                && optimState.pollcount <= nvars*2

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
                
                % Add vector to current point, fix to grid
                upollnew = bsxfun(@plus,u,vv);
                upollnew = periodCheck(upollnew,LB,UB,optimState);
                if options.ForcePollMesh
                    upollnew = force2grid(upollnew, optimState);
                end
                upollnew = uCheck(upollnew,optimState.TolMesh,optimState,0);
                
                % Add new poll points to polling set
                upoll = [upoll; upollnew];
                B = [B; Bnew];                    
            end
            
            % Cannot refill poll vector set, stop polling
            if isempty(upoll); break; end
            
            % Check whether it is time to refit the GP
            [refitgp_flag,unrelgp_flag,optimState] = IsRefitTime(optimState,options);
            if ~options.PollTraining && iter > 1; refitgp_flag = false; end
            
            % Rebuild GP local approximation if refitting or if at beginning of polling stage
            if refitgp_flag || optimState.pollcount == 0; gpstruct.post = []; end
            
            % Local GP approximation around polled points
            if isempty(gpstruct.post)
                [gpstruct,gptempflag] = gpTrainingSet(gpstruct, ...
                    options.gpMethod, ...
                    u, ...
                    upoll, ...
                    optimState, ...
                    options, ...
                    refitgp_flag);
                if refitgp_flag; refitted_flag = true; end
                gpexitflag = min(gptempflag,gpexitflag);
            end

            optimState = UpdateTarget(upollbest,fpollhyp,optimState,gpstruct,options);
                        
            % Evaluate acquisition function on poll vectors
            %--------------------------------------------------------------
            if options.AcqHedge
                [optimState.hedge,acqIndex,ymu,ys,fm,fs] = ...
                    acqPortfolio('acq',optimState.hedge,upoll,optimState.ftarget,0,gpstruct,optimState,options,SufficientImprovement);
                index = acqIndex(optimState.hedge.chosen);
            else
            %--------------------------------------------------------------
                % Batch evaluation of acquisition function on search set
                [z,~,ymu,ys,fm,fs] = feval(options.PollAcqFcn{:},upoll,optimState.ftarget,gpstruct,optimState,0);
                [~,index] = min(z);                
            end
            
            % Something went wrong, random vector
            if isempty(index) || isnan(index); index = randi(size(upoll,1)); end 
            
            % Compute probability that improvement at any location is 
            % less than SufficientImprovement (assumes independence --
            % conservative estimate towards continuing polling)
            gammaz = (optimState.ftarget - SufficientImprovement - fm)./fs;
            if numel(gpstruct.hyp) > 1
                gammaz = weightedsum(gpstruct.hypweight,gammaz,1);
            end
            if all(isfinite(gammaz)) && isreal(gammaz)
                fpi = 0.5*erfc(-gammaz/sqrt(2));
                fpi = sort(fpi,'descend');
                pless = prod(1-fpi(1:min(nvars,end)));
            else
                pless = 0;
                unrelgp_flag = 1;
            end
            
            % Consider whether to stop polling
            if ~options.CompletePoll            
                % Stop polling if last poll was good
                if goodpoll_flag
                    if unrelgp_flag
                        break;  % GP is unreliable, just stop polling                               
                    elseif pless > 1-options.TolPoI
                        break;  % Use GP prediction whether to stop polling
                    end
                else
                    % No good poll so far -- if GP is reliable, stop polling 
                    % if probability of improvement at any location is too low               
                    if ~unrelgp_flag && ... 
                            (options.ConsecutiveSkipping || lastskipped < iter-1) ...
                            && optimState.pollcount >= options.MinFailedPollSteps ...
                            && pless > 1-options.TolPoI
                        lastskipped = iter;
                        break;
                    end                    
                end
            end

            % Evaluate function and store value
            unew = upoll(index,:);
            [fpoll,optimState] = funlogger(funwrapper,unew,optimState,'iter');
            
            % Remove polled vector from set
            upoll(index,:) = [];

            % Save statistics of gp prediction
            optimState.gpstats = ...
                savegpstats(optimState.gpstats,fpoll,ymu(:,index),ys(:,index),gpstruct.hypweight);
            
            if optimState.UncertaintyHandling
                % Add just polled point to training set
                gpstruct = gpTrainingSet(gpstruct, ...
                    'add', ...
                    unew, ...
                    fpoll, ...
                    optimState, ...
                    options, ...
                    0);
                
                % Compute estimated function value at point
                [~,~,fpoll,fs2] = gppred(unew,gpstruct);
                if numel(gpstruct.hyp) > 1
                    fpoll = weightedsum(gpstruct.hypweight,fpoll,1);
                    fs2 = weightedsum(gpstruct.hypweight,fs2,1);
                end                
                fpollsd = sqrt(fs2);
            else
                fpollsd = 0;                
            end
            
            % Compute estimated improvement over incumbent
            PollImprovement = EvalImprovement(fval,fpoll,fsd,fpollsd,options.ImprovementQuantile);
            
            % Check if current point improves over best polled point so far
            if PollImprovement > PollBestImprovement 
                upollbest = unew;
                fpollbest = fpoll;
                fpollhyp = gpstruct.hyp;
                fpollbestsd = fpollsd;
                PollBestImprovement = PollImprovement;
                if PollBestImprovement > SufficientImprovement
                    goodpoll_flag = true;
                end
            end

            % Increase poll counter
            optimState.pollcount = optimState.pollcount + 1;
        end % Poll loop
        
        % Evaluate poll
        if (PollBestImprovement > 0 && options.SloppyImprovement) || ...
                PollBestImprovement > SufficientImprovement
            polldirection = find(abs(upollbest - ubest) > 1e-12,1); % The sign can be wrong for periodic variables (unused anyhow)            
            [ubest,fval,fsd,optimState,gpstruct] = UpdateIncumbent(ubest,fval,fsd,upollbest,fpollbest,fpollbestsd,optimState,gpstruct,options);
            u = ubest;
            pollmoved_flag = true;
        else
            pollmoved_flag = false;
        end

        if PollBestImprovement > SufficientImprovement
            % Successful poll, increase mesh size
            MeshSizeInteger = min(MeshSizeInteger + 1, options.MaxPollGridNumber);
            SuccessPoll_flag = true;
            optimState.usuccess = [optimState.usuccess; ubest];
            optimState.fsuccess = [optimState.fsuccess; fval];
        else
            % Failed poll, decrease mesh size
            MeshSizeInteger = MeshSizeInteger - 1;
            
            % Accelerated mesh reduction if stalling
            if options.AccelerateMesh && iter > options.AccelerateMeshSteps
                % Evaluate improvement in the last iterations
                HistoricImprovement = ...
                    EvalImprovement(optimState.iterList.fval(iter-options.AccelerateMeshSteps),fval,optimState.iterList.fsd(iter-options.AccelerateMeshSteps),fsd,options.ImprovementQuantile);
                if HistoricImprovement < options.TolFun
                    MeshSizeInteger = MeshSizeInteger - 1;
                end
            end            
            
            optimState.SearchSizeInteger = min(optimState.SearchSizeInteger, MeshSizeInteger*options.SearchGridMultiplier - options.SearchGridNumber);
            SuccessPoll_flag = false;
            
            % Profile plot of iteration
            if strcmpi(options.Plot,'profile') && ~isempty(gpstruct.x)
                % figure(iter);
                gpstruct.ftarget = optimState.ftarget;
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
        
        % Print iteration
        if prnt > 2
            if SuccessPoll_flag; string = 'Successful poll'; else string = 'Refine grid'; end
            action = [];
            if refitted_flag
                if isempty(action); action = 'Train'; else action = [action ', train']; end;
                if gpexitflag < 0; action = [action ' (failed)']; end
            end
            if lastskipped == iter; if isempty(action); action = 'Skip'; else action = [action ', skip']; end; end
            if optimState.UncertaintyHandling
                fprintf(displayFormat,iter,optimState.funccount,fval,fsd,MeshSize,string,action);                
            else
                fprintf(displayFormat,iter,optimState.funccount,fval,MeshSize,string,action);
            end
        end
                
        %H = fhess(@(xi_) gppred(xi_,gpstruct), gridunits(x,optimState), [], 'step', optimState.searchmeshsize)
        %H2 = fhess(@(x_) funwrapper(x_), x, [], 'step', optimState.searchmeshsize);        
        %[H; H2]
        % [inv(H2); optimState.Binv]
        
    end % Poll stage

    % Moved during the poll stage, need to retrain the GP
    if pollmoved_flag; gpstruct.post = []; end    
    
    %----------------------------------------------------------------------
    %% Finalize iteration
            
    % Scatter plot of iteration
    if strcmpi(options.Plot,'scatter')
        scatterplot(iter,ubest,fval,action,gpstruct,optimState,options);
    end
    
    fhyp = gpstruct.hyp;        % GP hyperparameters at end of iteration
        
    % Check termination conditions
    if optimState.funccount >= options.MaxFunEvals
        isFinished_flag = true;
        exitflag = 0;
        msg = 'Optimization terminated: reached maximum number of iterations OPTIONS.MaxIter.';
    end
    if iter >= options.MaxIter
        isFinished_flag = true;
        exitflag = 0;
        msg = 'Optimization terminated: reached maximum number of function evaluations OPTIONS.MaxFunEvals.';
    end
    if MeshSize < optimState.TolMesh
        isFinished_flag = true;
        exitflag = 1;
        msg = 'Optimization terminated: mesh size less than OPTIONS.TolMesh.';
    end
    if iter > options.TolStallIters
        HistoricImprovement = ...
            EvalImprovement(optimState.iterList.fval(iter-options.TolStallIters),fval,optimState.iterList.fsd(iter-options.TolStallIters),fsd,options.ImprovementQuantile);
        if HistoricImprovement < options.TolFun
            isFinished_flag = true;
            exitflag = 2;
            msg = 'Optimization terminated: change in the function value less than OPTIONS.TolFun.';
        end
    end
        
    % Store best points at the end of each iteration, or upon termination
    if DoPollStep_flag || isFinished_flag
        optimState.iterList.u(iter,:) = u;
        optimState.iterList.fval(iter,1) = fval;
        optimState.iterList.fsd(iter,1) = fsd;
        optimState.iterList.hyp{iter} = gpstruct.hyp;        
    end
    
    % Re-evaluate all points
    if DoPollStep_flag && optimState.UncertaintyHandling                        
        if iter > 1
            optimState = reevaluateIterList(optimState,gpstruct,options);
            
            % Update estimates of incumbent
            fval = optimState.iterList.fval(iter);
            fsd = optimState.iterList.fsd(iter);
            fhyp = optimState.iterList.hyp{iter};

            % Recompute improvement for all iterations
            ReImprovementList = EvalImprovement(fval,optimState.iterList.fval,fsd,optimState.iterList.fsd,options.ImprovementQuantile);
            [ReImprovement,index] = max(ReImprovementList);

            % Check if any point got better
            if ReImprovement > options.TolFun
                fval = optimState.iterList.fval(index);
                fsd = optimState.iterList.fsd(index);
                u = optimState.iterList.u(index,:);
                fhyp = optimState.iterList.hyp{index};
            end
        end
    end
    
    if isFinished_flag
        % Multiple starts (deprecated)
        if Restarts > 0 && optimState.funccount < options.MaxFunEvals
            isFinished_flag = false;
            MeshSizeInteger = 0;
            Restarts = Restarts - 1;
        end
    else
        % Iteration count is increased after the poll stage
        if DoPollStep_flag; iter = iter + 1; end        
    end

end

% Re-evaluate all best points (skip first iteration)
if optimState.UncertaintyHandling && iter > 1    
    optimState = reevaluateIterList(optimState,gpstruct,options);
        
    % Order by lowest probabilistic upper bound and choose best iterate
    SigmaMultiplier = sqrt(2).*erfcinv(2*options.FinalQuantile);    % Using inverted convention
    y = optimState.iterList.fval + SigmaMultiplier*optimState.iterList.fsd;
    [~,index] = min(y);

    % Best iterate
    fval = optimState.iterList.fval(index);
    fsd = optimState.iterList.fsd(index);
    u = optimState.iterList.u(index,:);
end

% Convert back to original space
x = origunits(u,optimState);

% Print final message
if prnt > 1
    fprintf('\n%s\n', msg);    
    if optimState.UncertaintyHandling
        fprintf('Estimated function value at minimum: %g ± %g (mean ± standard deviation).\n\n', fval, fsd);
    else
        fprintf('Function value at minimum: %g.\n\n', fval);
    end
end

if nargout > 3
    output.function = func2str(fun);    
    if optimState.UncertaintyHandling
        output.targettype = 'stochastic';
    else    
        output.targettype = 'deterministic';
    end    
    if all(isinf(LB)) && all(isinf(UB)) && isempty(nonbcon)
        output.problemtype = 'unconstrained';
    elseif isempty(nonbcon)
        output.problemtype = 'boundconstraints';
    else
        output.problemtype = 'nonboundconstraints';        
    end    
    output.iterations = iter;
    output.funccount = optimState.funccount;
    output.meshsize = optimState.meshsize;
    output.overhead = NaN;
    output.rngstate = rng;
    output.algorithm = 'Bayesian adaptive direct search';
    
    if ~isempty(optimState.nonbcon)
        output.maxconstraint = optimState.nonbcon(x);
    else
        output.maxconstraint = 0;
    end
    output.message = msg;
    
    % Return mean and SD of the estimated function value at the optimum
    output.fval = fval;
    output.fsd = fsd;
            
    % Return optimization struct (can be reused in future runs)
    if nargout > 4
        [~,optimState] = funlogger(funwrapper,u,optimState,'done');
    end
    
    % Compute total running time and fractional overhead
    optimState.totaltime = toc(t0);    
    output.overhead = optimState.totaltime / optimState.totalfunevaltime - 1;
end

end % Main BADS function

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
function [refitgp_flag,unrelgp_flag,optimState] = IsRefitTime(optimState,options)
%ISREFITTIME Check if time to refit the GP hyperparameters.

nvars = size(optimState.X,2);

% Check GP prediction statistics (track model reliability)
try
    unrelgp_flag = gppredcheck(optimState.gpstats,options.NormAlphaLevel);
catch
    unrelgp_flag = true;
end

% Refit more often early on (think if it can be improved)
if optimState.funccount < 200
    refitperiod = max(10, nvars*2);
else
    refitperiod = nvars*5;
end

refitgp_flag = optimState.lastfitgp < (optimState.funccount - options.MinRefitTime) && ...
    (optimState.gpstats.last >= refitperiod || unrelgp_flag) && ...
    optimState.funccount > nvars;

if refitgp_flag
    optimState.lastfitgp = optimState.funccount;
    % Reset GP prediction statistics
    Nsamples = max(1,options.gpSamples);
    optimState.gpstats = savegpstats([],[],[],[],ones(1,Nsamples)/Nsamples);
    unrelgp_flag = 0;
end

end

%--------------------------------------------------------------------------
function z = EvalImprovement(fbase,fnew,sbase,snew,q)
%EVALIMPROVEMENT Evaluate optimization improvement.
%   Z = EVALIMPROVEMENT(FBASE,FNEW) returns the improvement of FNEW over 
%   FBASE for a minimization problem (larger improvements are better).
%
%   Z = EVALIMPROVEMENT(FBASE,FNEW,SBASE,SNEW,Q) evaluates the quantile 
%   improvement for uncertain estimates of FBASE and FNEW with standard 
%   deviations, respectively, FNEW and FBASE. Q is the desired quantile.

if nargin < 3
    z = fbase - fnew;
else
    if q <= 0 || q >= 1
        error('Quantile Q for robust improvement should be greater than 0 and less than 1.');
    end
    mu = fbase - fnew;
    
    %----------------------------------------------------------------------
    % This needs to be corrected -- but for q=0.5 it does not matter
    sigma = sqrt(sbase.^2 + snew.^2);    
    x0 = -sqrt(2).*erfcinv(2*q);
    z = sigma.*x0 + mu;
    %----------------------------------------------------------------------
end

end

%--------------------------------------------------------------------------
function [unew,fvalnew,fsdnew,optimState,gpstruct] = UpdateIncumbent(uold,fvalold,fsdold,unew,fvalnew,fsdnew,optimState,gpstruct,options)
%UPDATEINCUMENT Move incumbent (current point) to a new point.

optimState.u = unew;
optimState.fval = fvalnew;
optimState.fsd = fsdnew;

% Update estimate of curvature (Hessian) -- not supported
if options.HessianUpdate && ~strcmpi(options.HessianMethod,'cmaes')
    updatehess(uold,fvalold,fsdold,unew,fvalnew,fsdnew,optimState,gpstruct,options);
end

end

%--------------------------------------------------------------------------
function optimState = UpdateTarget(ubest,hyp,optimState,gpstruct,options)
%UPDATETARGET Update optimization target (based on GP at incumbent).

if optimState.UncertaintyHandling || options.UncertainIncumbent
    gptemp = gpstruct;
    gptemp.hyp = hyp;
    [~,~,ftargetmu,ftargets2] = gppred(ubest,gptemp);
    ftargetmu = real(ftargetmu);
    if numel(gptemp.hyp) > 1
         ftargetmu = weightedsum(gptemp.hypweight,ftargetmu,1);
         ftargets2 = weightedsum(gptemp.hypweight,ftargets2,1);
    end    
    ftargets = sqrt(max(ftargets2,0));
    if ~isfinite(ftargetmu) || ~isreal(ftargets) || ~isfinite(ftargets)
        ftargetmu = optimState.fval;
        ftargets = optimState.fsd;
    end
    
    % Set optimization target slightly below current incumbent
    if ~options.AcqHedge || 1
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
%UPDATESEARCH Update search scale and statistics (for debugging purposes).

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

% Reset search factor at the end of each search
if optimState.searchcount == options.SearchNtry
    optimState.searchfactor = 1;        
end

end

%--------------------------------------------------------------------------
function optimState = reevaluateIterList(optimState,gpstruct,options)
%REEVALUATEITERLIST Re-evaluate function values at stored locations.

iter = optimState.iter;

% Return if no change since last re-evaluation
if optimState.lastreeval == optimState.funccount; return; end

% Loop over all iterations
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
    if numel(gpstruct.hyp) > 1
        fval = weightedsum(gpstruct.hypweight,fval,1);
        fs2 = weightedsum(gpstruct.hypweight,fs2,1);
    end    
    fsd = sqrt(fs2);

    optimState.iterList.fval(index) = fval;
    optimState.iterList.fsd(index) = fsd;
end

optimState.lastreeval = optimState.funccount;

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
function add2path()
%ADD2PATH Adds BADS subfolders to MATLAB path.

subfolders = {'acq','gpdef','gpml_fast','init','poll','search','utils','warp','gpml-matlab-v3.6-2015-07-07'};
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

end