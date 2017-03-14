function [u,fval,isFinished,optimState] = initmesh(u0,fval,funwrapper,SkipInitPoint,displayFormat,optimState,options)
%INITMESH Evaluate initial mesh

LB = optimState.LB;
UB = optimState.UB;
PLB = optimState.PLB;
PUB = optimState.PUB;
MeshSize = optimState.meshsize;

if SkipInitPoint && options.Ninit < 1
    error('If the initial point X0 is not specified, OPTIONS.Ninit needs to be > 0.');
end

% Only one function evaluation!
if ~SkipInitPoint && options.MaxFunEvals == 1
    u = u0;
    isFinished = 1;
    return;
end

% Additional initial points
if options.Ninit > 0

    if options.MaxFunEvals <= options.Ninit + double(~SkipInitPoint)
        %if strncmpi(options.Display,'all',3)
        %    fprintf('OPTIONS.MaxFunEvals <= OPTIONS.Ninit. Algorithm will only evaluate initial mesh.');
        %end
        Ninit = options.MaxFunEvals - double(~SkipInitPoint);
    else
        Ninit = options.Ninit;
    end

    % Call initialization function
    u1 = options.InitFcn(u0,LB,UB,PLB,PUB,Ninit,optimState,options);

    % Enforce periodicity
    u1 = periodCheck(u1,LB,UB,optimState);

    % Force points to be on the search grid
    u1 = force2grid(u1, optimState);

    % Ignore duplicates, vectors outside bounds or previously evaluated
    u1 = uCheck(u1,options.TolMesh,optimState,1);

    for i = 1:size(u1,1); [fval(i+1),optimState] = funlogger(funwrapper,u1(i,:),optimState,'iter'); end

    if all(~isfinite(fval)); error('Cannot find valid starting point.'); end

    [fval,idx] = min(fval);
    u1 = [u0; u1];
    u = u1(idx,:);

    if any(strcmpi(options.Display,{'iter','all'}))
        if options.UncertaintyHandling
            fprintf(displayFormat, 0, optimState.funccount, fval, NaN, MeshSize, 'Initial mesh', '');            
        else
            fprintf(displayFormat, 0, optimState.funccount, fval, MeshSize, 'Initial mesh', '');
        end
    end
else
    u = u0;
end

isFinished = 0;
