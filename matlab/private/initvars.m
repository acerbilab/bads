function [u0,LB,UB,PLB,PUB,MeshSizeInteger,MeshSize,TolMesh,optimState] = initvars(x0,LB,UB,PLB,PUB,optimState,options)
%INITVARS Initialize variables and transform coordinates

nvars = numel(x0);

if isscalar(LB); LB = LB*ones(1,nvars); end
if isscalar(UB); UB = UB*ones(1,nvars); end
if isscalar(PLB); PLB = PLB*ones(1,nvars); end
if isscalar(PUB); PUB = PUB*ones(1,nvars); end

if isempty(PLB) || isempty(PUB)
    warning('Plausible lower/upper bounds PLB and/or PUB not specified. Using hard upper/lower bounds instead.');
    if isempty(PLB); PLB = LB; end
    if isempty(PUB); PUB = UB; end
end

assert(all(isfinite(([PLB(:); PUB(:)]))), ...
    'Plausible interval bounds PLB and PUB need to be finite.');
    
% Grid parameters
MeshSizeInteger = 0;        % Mesh size in log base units
optimState.SearchSizeInteger = min(0,MeshSizeInteger*options.SearchGridMultiplier - options.SearchGridNumber);
MeshSize = options.PollMeshMultiplier^MeshSizeInteger;
optimState.meshsize = MeshSize;
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

u0 = force2grid(gridunits(x0(:)',optimState),optimState);
optimState.x0 = x0;             % Record starting point
optimState.u = u0;

% Put TOLMESH on space
TolMesh = options.PollMeshMultiplier^ceil(log(options.TolMesh)/log(options.PollMeshMultiplier));

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

% Setup covariance information
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

end