function options = setupoptions(nvars,defopts,options)
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
    'TolFun', 'TolNoise', 'Ninit', 'InitFcn', 'Restarts', 'CacheSize', 'FunValues', 'PeriodicVars', 'OutputFcn', ...
    'TolImprovement', 'ForcingExponent', 'PollMeshMultiplier', 'ForcePollMesh', 'IncumbentSigmaMultiplier', 'ImprovementQuantile', 'FinalQuantile', 'AlternativeIncumbent', 'AdaptiveIncumbentShift', 'FitnessShaping', 'WarpFunc', ...
    'NonlinearScaling', 'gpRescalePoll', 'PollMethod', 'Nbasis', ...
    'Nsearch', 'Nsearchiter', 'ESbeta', 'ESstart', 'SearchImproveFrac', 'SearchScaleSuccess', 'SearchScaleIncremental', 'SearchScaleFailure', 'SearchFactorMin', 'SearchMethod', ...
    'SearchGridNumber', 'MaxPollGridNumber', 'SearchGridMultiplier', 'SearchSizeLocked', 'SearchNtry', 'SearchMeshExpand', 'SearchMeshIncrement', 'SearchOptimize', ...
    'Ndata', 'MinNdata', 'BufferNdata', 'gpSamples', 'MinRefitTime', 'PollTraining', 'DoubleRefit', 'gpMeanPercentile', 'gpMeanRangeFun', ...
    'gpdefFcn', 'gpCluster','RotateGP','gpRadius','UseEffectiveRadius','PollAcqFcn', 'SearchAcqFcn', 'AcqHedge', 'CholAttempts', 'NoiseNudge', 'RemovePointsAfterTries', 'gpFixedMean', 'FitLik', 'gpSVGDiters', ...
    'UncertaintyHandling', 'NoiseObj', 'UncertainIncumbent', 'NoiseFinalSamples', 'TrustGPfinal', 'NoiseSize', 'MeshNoiseMultiplier', 'TolPoI', 'SkipPoll', 'ConsecutiveSkipping', 'SkipPollAfterSearch', 'CompletePoll', 'MinFailedPollSteps', 'NormAlphaLevel', ...
    'AccelerateMesh', 'AccelerateMeshSteps', 'SloppyImprovement', 'MeshOverflowsWarning', 'HessianUpdate', 'HessianAlternate', 'gpWarnings', ...
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
                error('bads:init', ...
                    'Cannot evaluate OPTIONS field "%s".', f{:});
            end
        end
    end
end

% Make cell arrays
cellfields = {'PollMethod','PollAcqFcn','SearchMethod','SearchAcqFcn'};
for f = cellfields
    if ischar(options.(f{:})) || isa(options.(f{:}), 'function_handle')
        options.(f{:}) = {options.(f{:})};
    end
end

% Check if MATLAB's Optimization Toolbox™ is available
if isempty(options.OptimToolbox)
    if exist('fmincon.m','file') && exist('fminunc.m','file') && exist('optimoptions.m','file') ...
            && license('test', 'optimization_toolbox')
        options.OptimToolbox = 1;
    else
        options.OptimToolbox = 0;
        warning('bads:noOptimToolbox', 'Could not find the Optimization Toolbox™. Using alternative optimization functions. This will slightly degrade performance. If you do not wish this message to appear, set OPTIONS.OptimToolbox = 0.');
    end
end

% Check options
if round(options.MaxFunEvals) ~= options.MaxFunEvals || options.MaxFunEvals <= 0
    error('OPTIONS.MaxFunEvals needs to be a positive integer.');
end

if options.ImprovementQuantile > 0.5
    warning('bads:excessImprovementQuantile', 'OPTIONS.ImprovementQuantile is greater than 0.5. This might produce unpredictable behavior. Set OPTIONS.ImprovementQuantile < 0.5 for conservative improvement.');
end

if ~isempty(options.NoiseSize) && options.NoiseSize(1) <= 0
    error('OPTIONS.NoiseSize, if specified, needs to be a positive scalar for numerical stability.');
end
    
