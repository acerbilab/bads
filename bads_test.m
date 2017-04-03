%% Test BADS

rng(0);
nvars = 6;                              % Number of dimensions
LB = -Inf(1,nvars);                     % Lower bound
UB = Inf(1,nvars);                      % Upper bound
PLB = -8*ones(1,nvars);                 % Plausible lower bound
PUB = 12*ones(1,nvars);                 % Plausible upper bound
x0 = (PUB-PLB).*rand(1,nvars) + PLB;    % Initial point

%% Ellipsoid (BPS)

display('Test with deterministic function (ellipsoid). Press any key to continue.');
pause;

fun = @(x) sum((x./(1:numel(x)).^2).^2);     % Objective function

options = bads('defaults');             % Default options
%options.UncertaintyHandling = 0;        % Deterministic function (determined at runtime if not specified)
%options.Plot = 'profile';               % Show profile during optimization
options.gpSamples = 1;
%options.gpdefFcn = '{@gpdefBads,''matern5'',1}';
%options.AcqHedge = 'on';
options.PollAcqFcn = '@acqThompson';
options.SearchAcqFcn = '@acqThompson';

rng(0);
[x,fval,exitflag,output,funValues,gpstruct] = bads(fun,x0,LB,UB,PLB,PUB,options);

display(['Final value: ' num2str(fval,'%.3f') ' (true value: 0.0), with ' num2str(output.funccount) ' fun evals.']);

%% Noisy sphere (BPS)

display('Test with noisy function (noisy sphere). Press any key to continue.');
pause;

fun = @(x) sum(x.^2) + randn();             % Noisy objective function

options = bads('defaults');             % Default options
%options.Plot = 'profile';              % Show profile during optimization
%options.UncertaintyHandling = 1;        % Activate noise handling (determined at runtime if not specified)
options.NoiseSize = 1;                  % Estimated noise magnitude
%options.gpSamples = 10;
%options.gpdefFcn = '{@gpdefBads,''matern5'',1}';

rng(0);
[x,fval,exitflag,output,funValues,gpstruct] = bads(fun,x0,LB,UB,PLB,PUB,options);

display(['Final value (not-noisy): ' num2str(sum(x.^2),'%.3f') ' (true value: 0.0) with ' num2str(output.funccount) ' fun evals.']);

%% Noisy sphere (FMINSEARCH)

display('Comparison with FMINSEARCH (noisy sphere). Press any key to continue.')

fun = @(x) sum(x.^2) + randn();             % Noisy objective function

options = optimset('Display','iter');

pause;

[x,fval,exitflag,output] = fminsearch(fun,x0,options);

display(['Final value (not-noisy): ' num2str(sum(x.^2),'%.3f') ' (true value: 0.0) with ' num2str(output.funcCount) ' fun evals.']);

%% Noisy sphere (GA)

if exist('ga.m','file')
    display('Comparison with GA (noisy sphere). Press any key to continue.')

    fun = @(x) sum(x.^2) + randn();             % Noisy objective function

    options = gaoptimset('Display','iter','Generations',19);
    
    pause;
    
    rng(0);
    [x,fval,exitFlag,output] = ga(fun,nvars,[],[],[],[],LB,UB,[],[],options);
    display(['Final value (not-noisy): ' num2str(sum(x.^2),'%.3f')  ' (true value: 0.0) with ' num2str(output.funccount) ' fun evals.']);
end