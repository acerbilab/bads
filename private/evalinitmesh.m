function [u,fval,isFinished,optimState,displayFormat] = evalinitmesh(u0,funwrapper,optimState,options,prnt)
%EVALINITMESH Evaluate initial mesh.

LB = optimState.LB;
UB = optimState.UB;
PLB = optimState.PLB;
PUB = optimState.PUB;
MeshSize = optimState.meshsize;

% If the objective fcn returns estimated noise, UncertaintyHandling is ON
if options.NoiseObj
    optimState.UncertaintyHandling = 1;
    if ~isempty(options.UncertaintyHandling) && ~options.UncertaintyHandling
        error('OPTIONS.NoiseObj is ON (objective fcn with estimated noise output) but OPTIONS.UncertaintyHandling is OFF.');
    end
else
    optimState.UncertaintyHandling = options.UncertaintyHandling;
end

fval = [];

% Evaluate non-bound constraints
if ~isempty(optimState.nonbcon)
    c = optimState.nonbcon(optimState.x0);  % Evaluate constraints
    if c > 0; fval = NaN; end
end

% Evaluate starting point (if specified)
if isempty(fval)
    [fval,optimState] = funlogger(funwrapper,u0,optimState,'iter');
end
optimState.fval = fval;

% If UncertaintyHandling is not specified, test if function is noisy
if isempty(optimState.UncertaintyHandling)
    uncertaintyTest = 'Uncertainty test';
    x0 = origunits(u0,optimState);
    fval2 = funwrapper(x0);
    optimState.funccount = optimState.funccount + 1;    % +1 fcn eval to count
    if abs(fval - fval2) > options.TolNoise
        optimState.UncertaintyHandling = 1;
        if prnt > 0
            fprintf('Beginning optimization of a STOCHASTIC objective fcn.\n');
        end
    else
        optimState.UncertaintyHandling = 0;
        if prnt > 1
            fprintf('Beginning optimization of a DETERMINISTIC objective fcn.\n');
        end
    end
else
    uncertaintyTest = '';    
end

% Display format depends whether the objective is noisy
if prnt > 2
    if optimState.UncertaintyHandling
        displayFormat = ' %5.0f       %5.0f    %12.6g    %12.6g    %12.6g    %20s      %s\n';
        fprintf(' Iteration    f-count      E[f(x)]        SD[f(x)]           MeshScale          Method          Actions\n');
    else
        displayFormat = ' %5.0f       %5.0f    %12.6g    %12.6g    %20s      %s\n';
        fprintf(' Iteration    f-count         f(x)           MeshScale          Method          Actions\n');
    end
else
    displayFormat = [];
end

if prnt > 2
    if optimState.UncertaintyHandling
        fprintf(displayFormat, 0, optimState.funccount, fval, NaN, MeshSize, '', uncertaintyTest);
    else
        fprintf(displayFormat, 0, optimState.funccount, fval, MeshSize, '', uncertaintyTest);
    end
end

% Only one function evaluation!
if options.MaxFunEvals == 1
    u = u0;
    isFinished = 1;
    return;
end

% If dealing with a noisy function, use a large initial mesh
if optimState.UncertaintyHandling
    options.Ninit = min(max(20,options.Ninit),options.MaxFunEvals);
end

% Additional initial points
if options.Ninit > 0

    % Evaluate initial points but not more than OPTIONS.MaxFunEvals
    Ninit = min(options.Ninit, options.MaxFunEvals - 1);

    % Call initialization function
    u1 = options.InitFcn(u0,LB,UB,PLB,PUB,Ninit,optimState,options);

    % Enforce periodicity
    u1 = periodCheck(u1,LB,UB,optimState);

    % Force points to be on the search grid
    u1 = force2grid(u1, optimState);

    % Ignore duplicates, vectors outside bounds or previously evaluated
    u1 = uCheck(u1,optimState.TolMesh,optimState,1);

    % Evaluate all points
    for i = 1:size(u1,1)
        [fval(i+1),optimState] = funlogger(funwrapper,u1(i,:),optimState,'iter');
    end

    [fval,idx] = min(fval);
    u1 = [u0; u1];
    u = u1(idx,:);

    if prnt > 2
        if optimState.UncertaintyHandling
            fprintf(displayFormat, 0, optimState.funccount, fval, NaN, MeshSize, 'Initial mesh', '');            
        else
            fprintf(displayFormat, 0, optimState.funccount, fval, MeshSize, 'Initial mesh', '');
        end
    end
else
    u = u0;
end

optimState.fval = fval;
isFinished = 0;
