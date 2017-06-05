%BADS_EXAMPLES Examples for Bayesian Adaptive Direct Search
%
%  Example 1: Basic usage
%  Example 2: Non-bound constraints
%  Example 3: Noisy objective function
%  Example 4: Extra noisy objective function
%  Example 5: Periodic function
%  Example 6: Extended usage
%
%  Note: after installation, run 
%    bads('test') 
%  to check that everything is working correctly.
%
%  For any question, check out the FAQ: 
%  https://github.com/lacerbi/bads/wiki
%
%  See also BADS.

% Luigi Acerbi 2017

display('Running a number of examples usage for Bayesian Adaptive Direct Search (BADS).');
display('Open ''bads_examples.m'' to see additional comments and instructions.');


%% Example 1: Basic usage

% Simple usage of BADS on Rosenbrock's banana function in 2D
% (see https://en.wikipedia.org/wiki/Rosenbrock_function).
% 
% We specify wide hard bounds and tighter plausible bounds that (hopefully) 
% contain the solution. Plausible bounds represent your best guess at 
% bounding the region where the solution might lie.

x0 = [0 0];                 % Starting point
lb = [-20 -20];             % Lower bounds
ub = [20 20];               % Upper bounds
plb = [-5 -5];              % Plausible lower bounds
pub = [5 5];                % Plausible upper bounds

% Screen display
fprintf('\n');
display('*** Example 1: Basic usage');
display('  Simple usage of BADS on <a href="https://en.wikipedia.org/wiki/Rosenbrock_function">Rosenbrock''s banana function</a> in 2D.');
display('  Press any key to continue.'); fprintf('\n');
pause;

% Run BADS, which returns the minimum X and its value FVAL.
[x,fval] = bads(@rosenbrocks,x0,lb,ub,plb,pub)

display('The true global minimum is at X = [1,1], where FVAL = 0.');

% Note that BADS by default does not aim for extreme numerical precision 
% (e.g., beyond the 2nd or 3rd decimal place), since in realistic 
% model-fitting problems such a resolution is typically pointless.


%% Example 2: Non-bound constraints

% We test BADS by forcing the solution to be within a circle with radius 1.
% Since we know the optimization region, we set tight hard bounds around it
% to further help the search.

x0 = [0,0];                 % Starting point
lb = [-1 -1];   ub = [1 1]; % Hard bounds only

% Note that BADS will complain because the plausible bounds are not 
% specified explicitly (it will use LB and UB instead). Generally, you want
% to specify both hard and plausible bounds.

% Non-bound constraints are violated outside the unit circle
nonbcon = @(x) sum(x.^2,2) > 1;

% Note that NONBCON requires a matrix input. Suppose we want to write the 
% above case without using SUM. We would have:
% nonbcon = @(x) (x(:,1).^2 + x(:,2).^2) > 1;   % Correct
% nonbcon = @(x) (x(1).^2 + x(2).^2) > 1;       % Wrong! not matrix input

% Screen display
fprintf('\n');
display('*** Example 2: Non-bound constraints');
display('  As before, but we force the input to stay in a circle with unit radius.');
display('  BADS will complain because the plausible bounds are not specified explicitly.');
display('  Press any key to continue.'); fprintf('\n');
pause;

% Run BADS with both bound and non-bound constraints
[x,fval] = bads(@rosenbrocks,x0,lb,ub,[],[],nonbcon)

% Alternatively, the following instructions would make BADS happier 
% (in a realistic model-fitting scenario, we recommend, whenever possible, 
% to specify plausible bounds which are tighter than the hard bounds).
%
% plb = lb; pub = ub;
% [x,fval] = bads(@rosenbrocks,x0,lb,ub,plb,pub,nonbcon)

display('The true global minimum under these constraints is at X = [0.786,0.618], where FVAL = 0.046.');


%% Example 3: Noisy objective function

% We test BADS on a noisy function. For the purpose of this test, we
% manually add unit Gaussian noise to the 'sphere' (quadratic) function.

% Define the noisy objective as a function handle
noisyfun = @(x) sum(x.^2,2) + randn(size(x,1),1);

x0 = [-3 -3];       % For a change, we start farther away from the solution
% This time to help the search we set tighter bounds
lb = [-5 -5];   ub = [5 5];
plb = [-2 -2];  pub = [2 2];

% Screen display
fprintf('\n');
display('*** Example 3: Noisy objective function');
display('  We test BADS on a noisy quadratic function with unit Gaussian noise.');
display('  Press any key to continue.'); fprintf('\n');
pause;

% Run BADS on the noisy function
[x,fval,exitflag,output] = bads(noisyfun,x0,lb,ub,plb,pub);
x
fval
display(['The true, noiseless value of the function at X is ' num2str(sum(x.^2,2)) '.']);
display('The true global minimum is at X = [0,0], where FVAL = 0.');

% FVAL in this case is an *estimate* of the function at X, obtained by
% averaging ten function evaluations. These values can be found in the
% OUTPUT structure, together with additional information about the optimization.

display('The returned OUTPUT structure is:');
output

% Note that the fractional overhead of BADS reported in OUTPUT is astronomical.
% The reason is that the objective function we are using is analytical and 
% extremely fast, which is not what BADS is designed for. 
% In a realistic scenario, the objective function will be moderately costly
% (e.g., more than 0.1 s per function evaluation), and the fractional 
% overhead should be less than 1.


%% Example 4: Extra noisy objective function

% We test BADS on a particularly noisy function and look at some options.

% Define noisy objective with substantial input-dependent noise
noisyfun = @(x) sum(x.^2,2) + (3 + 0.1*sqrt(sum(x.^2,2))).*randn(size(x,1),1);

% For this optimization, we explicitly tell BADS that the objective is
% noisy (it is not necessary, but it is a good habit); and also specify a 
% rough estimate for the value of the noise in a neighborhood of the solution.
% Finally, we tell BADS to use more samples to estimate FVAL at the end.

options = [];                       % Reset the OPTIONS struct
options.UncertaintyHandling = 1;    % Tell BADS that the objective is noisy
options.NoiseSize           = 5;    % Estimate of noise
options.NoiseFinalSamples   = 100;  % # samples to estimate FVAL at the end 
                                    % (default would be 10)
x0 = [-3 -3];
lb = [-5 -5];   ub = [5 5];
plb = [-2 -2];  pub = [2 2];

% Screen display
fprintf('\n');
display('*** Example 4: Extra noisy function');
display('  We test BADS on a particularly noisy function.');
display('  Press any key to continue.'); fprintf('\n');
pause;

% Run BADS on the noisy function
[x,fval,exitflag,output] = bads(noisyfun,x0,lb,ub,plb,pub,[],options);
x
fval
display(['The true, noiseless value of the function at X is ' num2str(sum(x.^2,2)) '.']);
display('The true global minimum is at X = [1,1], where FVAL = 0.');
display('Due to the elevated level of noise, we do not necessarily expect high precision in the solution.');


%% Example 5: Objective function with periodic dimensions

% We test BADS on a function with a subset of periodic dimensions.

% This function is periodic along the third and fourth dimension, with
% periods respectively 4 and 2.
periodicfun = @(x) rosenbrocks(x(:,1:2)) + cos(x(:,3)*pi/2) + cos(x(:,4)*pi) + 2;

x0 = [-3 -3 -1 -1];
% We specify the periodic bounds via hard bounds
lb = [-10 -5 -2 -1];
ub = [5 10 2 1];

plb = [-2 -2 -2 -1];
pub = [2 2 2 1];

options = [];                       % Reset the OPTIONS struct
options.PeriodicVars = [3 4];       % The 3rd and 4th variables are periodic

% Screen display
fprintf('\n');
display('*** Example 5: Objective function with periodic dimensions');
display('  We test BADS on a function with some periodic inputs.');
display('  Press any key to continue.'); fprintf('\n');
pause;

[x,fval,exitflag,output] = bads(periodicfun,x0,lb,ub,plb,pub,[],options);
x
fval
display('The true global minimum is at X = [1,1,±2,±1], where FVAL = 0.');


%% Example 6: Extended usage

% Extended usage of BADS that shows some additional options.

% Function handle for function with multiple input arguments (e.g., here
% we add a translation of the input; but more in general you could be 
% passing additional data to your objective function).
fun = @(x,mu) rosenbrocks(bsxfun(@plus, x, mu)); 

% This will translate the Rosenbrock fcn such that the global minimum is at zero
mu = [1 1 1 1];

% We now set bounds using also fixed variables
% (2nd and 4th variable are fixed by setting all bounds and X0 equal)

plb = [-2 0 -2 0];             % Plausible lower bounds
pub = [2 0 2 0];               % Plausible upper bounds
lb = [-20 0 -5 0];             % Hard lower bounds
ub = [20 0 5 0];               % Hard upper bounds

% Random starting point inside plausible region. In a typical optimization
% scenario, you will repeat the optimization from different starting
% points (ideally 10 or more), possibly drawn this way, and take the best
% result.
x0 = plb + (pub-plb).*rand(1,numel(plb));

options = bads('defaults');             % Get a default OPTIONS struct
options.MaxFunEvals         = 50;       % Very low budget of function evaluations
options.Display             = 'final';   % Print only basic output ('off' turns off)
options.UncertaintyHandling = 0;        % We tell BADS that the objective is deterministic

% Custom output function (return FALSE to continue, TRUE to stop optimization)
options.OutputFcn           = @(x,optimState,state) ~isfinite(fprintf('%s %d... ', state, optimState.iter));

% Screen display
fprintf('\n');
display('*** Example 6: Extended usage');
display('  Extended usage of BADS with additional options and no detailed display.');
display('  Press any key to continue.'); fprintf('\n');
pause;

% Run BADS, passing MU as additional (fixed) input argument for FUN
[x,fval,exitflag,output] = bads(fun,x0,lb,ub,plb,pub,[],options,mu);

% The following line of code would do the same using an anonymous function
% [x,fval,exitflag] = bads(@(x) fun(x,mu),x0,lb,ub,plb,pub,[],options);

x
fval
display('The true global minimum is at X = [0,0,0,0], where FVAL = 0.');
exitflag
display('EXITFLAG of 0 means that the maximum number of function evaluations has been reached.');
fprintf('\n');
display('For this optimization we used the following OPTIONS:')
options
display('Type ''help bads'' for additional documentation on BADS, or consult the <a href="https://github.com/lacerbi/bads">Github page</a> or <a href="https://github.com/lacerbi/bads/wiki">online FAQ</a>.');