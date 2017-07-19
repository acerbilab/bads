function [u0,LB,UB,PLB,PUB,MeshSizeInteger,optimState] = setupvars(x0,LB,UB,PLB,PUB,optimState,nonbcon,options,prnt)
%INITVARS Initialize variables and transform coordinates.

nvars = numel(x0);

% Test bound order
if ~all(LB <= PLB & PLB < PUB & PUB <= UB)
    error('Bound vectors do not respect the order: LB <= PLB < PUB <= UB.');
end

% Check NONBCON
if ~isempty(nonbcon)
    conerr = false;
    try
        y = nonbcon([PLB; PUB]);
        if size(y,1) ~= 2 || size(y,2) ~= 1
            conerr = true;
        end
    catch
        conerr = true;        
    end
    if conerr
        error('NONBCON should be a function handle that takes a matrix X as input and returns a column vector of bound violations.');        
    end
end

% Gentle warning for infinite bounds
nInfs = sum(isinf([LB(:);UB(:)]));
if nInfs > 0
    if prnt > 0
        if nInfs == 2*nvars
            fprintf('Caution: Detected fully unconstrained optimization.\n');
        else
            fprintf('Caution: Detected %d infinite bound(s).\n', nInfs);
        end
        fprintf('Infinite lower/upper bounds are deprecated. Support might be removed in future versions of BADS.\n\n');
    end
end

% Grid parameters
MeshSizeInteger = 0;        % Mesh size in log base units
optimState.SearchSizeInteger = min(0,MeshSizeInteger*options.SearchGridMultiplier - options.SearchGridNumber);
optimState.meshsize = options.PollMeshMultiplier^MeshSizeInteger;
optimState.searchmeshsize = options.PollMeshMultiplier.^optimState.SearchSizeInteger;
optimState.scale = 1;

% Compute transformation of variables

if options.NonlinearScaling
    % Let TRANSVARS decide on nonlinear rescaling of variables
    logflag = NaN(1,nvars);
    if ~isempty(options.PeriodicVars)   % Never transform periodic variables
        logflag(options.PeriodicVars) = 0;
    end
else
    logflag = zeros(1,nvars);           % Do not allow nonlinear rescaling
end
optimState.trinfo = transvars(nvars,LB,UB,PLB,PUB,logflag); 

% Convert bounds to standardized space
LB = optimState.trinfo.lb;
UB = optimState.trinfo.ub;
PLB = optimState.trinfo.plb;
PUB = optimState.trinfo.pub;
optimState.LB = LB;
optimState.UB = UB;
optimState.PLB = PLB;
optimState.PUB = PUB;

% Bounds for search mesh
optimState.LBsearch = force2grid(LB,optimState);
optimState.LBsearch(optimState.LBsearch < optimState.LB) = ...
    optimState.LBsearch(optimState.LBsearch < optimState.LB) + optimState.searchmeshsize;
optimState.UBsearch = force2grid(UB,optimState);
optimState.UBsearch(optimState.UBsearch > optimState.UB) = ...
    optimState.UBsearch(optimState.UBsearch > optimState.UB) - optimState.searchmeshsize;

% Starting point in grid coordinates
if any(~isfinite(x0))   % Invalid/not provided starting point
    if prnt > 0
        fprintf('Initial starting point is invalid or not provided. Starting from center of plausible region.\n');
    end
    u0 = force2grid(0.5*(PLB + PUB),optimState);    % Midpoint
    x0 = origunits(u0,optimState);
else
    u0 = force2grid(gridunits(x0,optimState),optimState);
end

% Adjust points that fall outside bounds due to gridization
u0(u0 < optimState.LB) = u0(u0 < optimState.LB) + optimState.searchmeshsize;
u0(u0 > optimState.UB) = u0(u0 > optimState.UB) - optimState.searchmeshsize;

% Test that starting point X0 (U0 in transformed space) is within bounds
LB_orig = optimState.trinfo.oldbounds.lb;
UB_orig = optimState.trinfo.oldbounds.ub;
if ~all(x0 <= UB_orig & x0 >= LB_orig & u0 <= UB & u0 >= LB)
    error('Initial starting point is not within the hard bounds LB and UB.');
end

optimState.x0 = x0; % Record starting point (original coordinates)
optimState.u = u0;

% Put TOLMESH on space
optimState.TolMesh = options.PollMeshMultiplier^ceil(log(options.TolMesh)/log(options.PollMeshMultiplier));

% Periodic variables
optimState.periodicvars = false(1, nvars);
if ~isempty(options.PeriodicVars)
    optimState.periodicvars(options.PeriodicVars) = true;
    for d = find(optimState.periodicvars)
        if ~isfinite(LB(d)) || ~isfinite(UB(d))
            error('Periodic variables need to have finite lower and upper bounds.');
        end
    end
end

% Report variable transformation
if any(optimState.trinfo.logct) && prnt > 0
    fprintf('Variables (index) internally transformed to log coordinates: %s.\n',mat2str(find(optimState.trinfo.logct)));
end
if any(optimState.periodicvars) && prnt > 0
    fprintf('Variables (index) defined with periodic boundaries: %s.\n',mat2str(find(optimState.periodicvars)));
end


% Setup covariance information (unused)
if options.HessianUpdate
    if strcmpi(options.HessianMethod,'cmaes')
        D = ((PUB-PLB)./optimState.scale).^2/12;
        optimState.Binv = diag(D/sum(D));
        optimState.C = sqrt(optimState.Binv);
    else
        optimState.Binv = eye(nvars);        
        optimState.C = eye(nvars);    
    end
    optimState.grad = ones(nvars,1);
    optimState.violations = 0;
    optimState.B = inv(optimState.Binv);
end

% Import prior function evaluations
if ~isempty(options.FunValues)
    if ~isfield(options.FunValues,'X') || ~isfield(options.FunValues,'Y')
        error('bads:funValues', ...
            'The ''FunValues'' field in OPTIONS needs to have a X and a Y field (respectively, inputs and their function values).');
    end
        
    X = options.FunValues.X;
    Y = options.FunValues.Y;
    if size(X,1) ~= size(Y,1)
        error('X and Y arrays in the OPTIONS.FunValues need to have the same number of rows (each row is a tested point).');        
    end
    
    if ~all(isfinite(X(:))) || ~all(isfinite(Y(:))) || ~isreal(X) || ~isreal(Y)
        error('X and Y arrays need to be finite and real-valued.');                
    end    
    if ~isempty(X) && size(X,2) ~= nvars
        error('X should be a matrix of tested points with the same dimensionality as X0 (one input point per row).');
    end
    if ~isempty(Y) && size(Y,2) ~= 1
        error('Y should be a vertical array of function values (one function value per row).');
    end
    
    optimState.X = X;
    optimState.Y = Y;    
    
    % Heteroskedastic noise
    if isfield(options.FunValues,'S')
        S = options.FunValues.S;
        if size(S,1) ~= size(Y,1)
            error('X, Y, and S arrays in the OPTIONS.FunValues need to have the same number of rows (each row is a tested point).');        
        end    
        if ~all(isfinite(S)) || ~isreal(S) || ~all(S > 0)
            error('S array needs to be finite, real-valued, and positive.');
        end
        if ~isempty(S) && size(S,2) ~= 1
            error('S should be a vertical array of estimated function SD values (one SD per row).');
        end
        optimState.S = S;        
    end    
    
end

%% Initialize OPTIMSTATE variables

optimState.searchfactor =   1;
optimState.sdlevel      = options.IncumbentSigmaMultiplier;
optimState.searchcount  = options.SearchNtry;       % Skip search at first iteration
optimState.lastreeval   = -Inf;                     % Last time function values were re-evaluated
optimState.lastfitgp    = -Inf;                     % Last fcn evaluation for which the gp was trained
optimState.meshoverflows = 0;                       % Number of attempted mesh expansions when already at maximum size

% List of points at the end of each iteration
optimState.iterList.u = [];
optimState.iterList.fval = [];
optimState.iterList.fsd = [];
optimState.iterList.fhyp = [];

% Create vector of ES weights (only for searchES)
es_iter = options.Nsearchiter;
es_mu = options.Nsearch/es_iter;
es_lambda = es_mu;
optimState.es = ESupdate(es_mu,es_lambda,es_iter);

% Hedge struct
optimState.hedge = [];

