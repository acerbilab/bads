function output = bads_output(fun,optimState,options,LB,UB,nonbcon,iter,x,msg,yval_vec,ysd_vec,fval,fsd,totaltime,bads_version)
%BADS_OUTPUT Create OUTPUT struct for BADS.

output.function = func2str(fun);    
if optimState.UncertaintyHandling
    if options.SpecifyTargetNoise
        output.targettype = 'stochastic (specified noise)';
    else
        output.targettype = 'stochastic';
    end
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
output.version = bads_version;

if ~isempty(optimState.nonbcon)
    output.maxconstraint = optimState.nonbcon(x);
else
    output.maxconstraint = 0;
end
output.message = msg;

% Observed function value(s) at optimum (possibly multiple samples)
output.yval = yval_vec;
if ~isempty(ysd_vec)
    output.ysd = ysd_vec;
end

% Return mean and SD of the estimated function value at the optimum
output.fval = fval;
output.fsd = fsd;

% Compute total running time and fractional overhead
optimState.totaltime = totaltime;    
output.overhead = optimState.totaltime / optimState.totalfunevaltime - 1;

end