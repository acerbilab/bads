function [x,fval,exitflag,output,optimState,gpstruct] = fixedbads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options,fixidx,nout)
%FIXEDBADS Run BADS with reduced dimensionality due to fixed variables.

nvars = size(x0,2);

% Initialize outputs
fval = [];
exitflag = [];
output = [];
optimState = [];
gpstruct = [];

freeidx = ~fixidx;
fixedvars = LB(fixidx);

% Update periodic variables, if present in OPTIONS
if isfield(options,'PeriodicVars') && ~isempty(eval(options.PeriodicVars))
    pervars = false(1,nvars);
    pervars(eval(options.PeriodicVars)) = true;
    pervars(fixidx) = [];
    options.PeriodicVars = find(pervars);
end

% Fix input FUNVALUES, if present
if isfield(options,'FunValues') && ~isempty(options.FunValues) && ...
        isfield(options.FunValues,'X') && ~isempty(options.FunValues.X)
    options.FunValues.X = options.FunValues.X(:,freeidx);
end

% Downlift
x0 = x0(freeidx);
LB = LB(freeidx);
UB = UB(freeidx);
PLB = PLB(freeidx);
PUB = PUB(freeidx);

switch nout
    case 1
        x_free = bads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options);    
    case 2
        [x_free,fval] = bads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options);
    case 3
        [x_free,fval,exitflag] = bads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options);
    case 4
        [x_free,fval,exitflag,output] = bads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options);
    case 5
        [x_free,fval,exitflag,output,optimState] = bads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options);
    case 6
        [x_free,fval,exitflag,output,optimState,gpstruct] = bads(fun_fix,x0,LB,UB,PLB,PUB,nonbcon_fix,options);
end

% Uplift solution
x = zeros(1,nvars);
x(freeidx) = x_free;
x(fixidx) = fixedvars;

if nout > 4
    % Uplift points in OPTIMSTATE
    nx = size(optimState.X,1);
    Xtemp = zeros(nx,nvars);
    Xtemp(:,freeidx) = optimState.X;
    Xtemp(:,fixidx) = repmat(fixedvars, [nx 1]);
    optimState.X = Xtemp;
end

end