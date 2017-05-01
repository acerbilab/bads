function [u0,LB,UB,PLB,PUB,MeshSizeInteger,optimState] = setupvars(x0,LB,UB,PLB,PUB,optimState,options,prnt)
%INITVARS Initialize variables and transform coordinates.

nvars = numel(x0);

% Expand scalar inputs to full vectors
if isscalar(LB); LB = LB*ones(1,nvars); end
if isscalar(UB); UB = UB*ones(1,nvars); end
if isscalar(PLB); PLB = PLB*ones(1,nvars); end
if isscalar(PUB); PUB = PUB*ones(1,nvars); end

if isempty(PLB) || isempty(PUB)
    warning('Plausible lower/upper bounds PLB and/or PUB not specified. Using hard upper/lower bounds instead.');
    if isempty(PLB); PLB = LB; end
    if isempty(PUB); PUB = UB; end
end

% Test that all vectors have the same length
if any([numel(LB),numel(UB),numel(PLB),numel(PUB)] ~= nvars)
    error('All input vectors (X0, LB, UB, PLB, PUB), if specified, need to have the same size.');
end

% Test that all vectors are row vectors
if ~isvector(x0) || ~isvector(LB) || ~isvector(UB) || ~isvector(PLB) || ~isvector(PUB) ...
        || size(x0,1) ~= 1 || size(LB,1) ~= 1 || size(UB,1) ~= 1 || size(PLB,1) ~= 1 || size(PUB,1) ~= 1
    error('All input vectors should be row vectors.');
end

% Test that plausible bounds are finite
if ~all(isfinite(([PLB, PUB]))) 
    error('Plausible interval bounds PLB and PUB need to be finite.');
end

% Test that all vectors are real-valued
if ~isreal([x0; LB; UB; PLB; PUB])
    error('All input vectors should be real-valued.');
end

% Test bound order
if ~all(LB <= PLB & PLB < PUB & PUB <= UB)
    error('Bound vectors do not respect the order: LB <= PLB < PUB <= UB.');
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
    if ~all(isfinite(X)) || ~all(isfinite(Y)) || ~isreal(X) || ~isreal(Y)
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

