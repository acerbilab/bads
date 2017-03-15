function options = initoptions(nvars,defopts,options)
%INITOPTIONS Initialize OPTIONS struct.

% Assign default values to OPTIONS struct
for f = fieldnames(defopts)'
    if ~isfield(options,f{:}) || isempty(options.(f{:}))
        options.(f{:}) = defopts.(f{:});
    end
end

% Remove comments and trailing empty spaces from options fields
for f = fieldnames(options)'
    if ischar(options.(f{:}))
        idx = find(options.(f{:}) == '%',1);
        if ~isempty(idx); options.(f{:})(idx:end) = []; end        
        idx = find(options.(f{:}) ~= ' ',1,'last');
        if ~isempty(idx); options.(f{:})(idx+1:end) = []; end                
    end
end

% OPTIONS fields that need to be evaluated
evalfields = {'Debug', 'MaxIter', 'MaxFunEvals', 'TolMesh', 'TolStallIters', ...
    'TolFun', 'Ninit', 'InitFcn', 'Restarts', 'CacheSize', 'FunValues', 'PeriodicVars' ...
    'TolImprovement', 'ForcingExponent', 'PollMeshMultiplier', 'IncumbentSigmaMultiplier', 'AlternativeIncumbent', 'AdaptiveIncumbentShift', 'FitnessShaping', 'WarpFunc', ...
    'NonlinearScaling', 'gpRescalePoll', 'PollMethod', 'Nbasis', ...
    'Nsearch', 'Nsearchiter', 'ESbeta', 'ESstart', 'SearchImproveFrac', 'SearchScaleSuccess', 'SearchScaleIncremental', 'SearchScaleFailure', 'SearchFactorMin', 'SearchMethod', ...
    'SearchGridNumber', 'MaxPollGridNumber', 'SearchGridMultiplier', 'SearchSizeLocked', 'SearchNtry', 'SearchMeshExpand', 'SearchMeshIncrement', 'SearchOptimize', ...
    'Ndata', 'MinNdata', 'BufferNdata', 'gpSamples', 'gpMarginalize', 'MinRefitTime', 'PollTraining', 'DoubleRefit', 'gpMeanPercentile', 'gpMeanRangeFun', ...
    'gpdefFcn', 'gpCluster','RotateGP','gpRadius','UseEffectiveRadius','PollAcqFcn', 'SearchAcqFcn', 'AcqHedge', 'CholAttempts', 'NoiseNudge', 'RemovePointsAfterTries', 'gpFixedMean', 'FitLik', ...
    'UncertaintyHandling', 'NoiseObj', 'UncertainIncumbent', 'NoiseSize', 'MeshNoiseMultiplier', 'TolPoI', 'SkipPoll', 'ConsecutiveSkipping', 'SkipPollAfterSearch', 'MinFailedPollSteps', 'NormAlphaLevel', ...
    'AccelerateMesh', 'AccelerateMeshSteps', 'SloppyImprovement', 'HessianUpdate', 'HessianAlternate', ...
    'HedgeGamma','HedgeBeta','HedgeDecay', ...
    'TrueMinX', 'OptimToolbox' ...
    };

% Evaluate string options
for f = evalfields
    if ischar(options.(f{:}))
        try
            options.(f{:}) = eval(options.(f{:}));
        catch
            try
                options.(f{:}) = evalbool(options.(f{:}));
            catch
                error('bps:init', ...
                    'Cannot evaluate OPTIONS field "%s".', f{:});
            end
        end
    end
end

% Make cell arrays
cellfields = {'PollMethod','PollAcqFcn','SearchMethod','SearchAcqFcn'};
for f = cellfields
    if ischar(options.(f{:})) || isa(options.(f{:}), 'function_handle');
        options.(f{:}) = {options.(f{:})};
    end
end

% Change options for uncertainty handling
if options.UncertaintyHandling
    options.TolStallIters = 2*options.TolStallIters;
    options.Ndata = max(200,options.Ndata);
    options.MinNdata = 2*options.MinNdata;
    options.Ninit = min(max(20,options.Ninit),options.MaxFunEvals);
    %options.gpMeanPercentile = 50;
    options.MinFailedPollSteps = Inf;
    options.MeshNoiseMultiplier = 0;
end

% Check if MATLAB's Optimization Toolbox™ is available
if isempty(options.OptimToolbox)
    if exist('fmincon.m','file') && exist('fminunc.m','file') && exist('optimoptions.m','file') ...
            && license('test', 'optimization_toolbox')
        options.OptimToolbox = 1;
    else
        options.OptimToolbox = 0;
        warning('Could not find the Optimization Toolbox™. Using alternative optimization functions. This will slightly degrade performance. If you do not wish this message to appear, set OPTIONS.OptimToolbox = 0.');
    end    
end

end