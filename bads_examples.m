%BADS_EXAMPLES Examples for Bayesian Adaptive Direct Search
%
%  Example 1: Basic usage
%  Example 2: Non-bound constraints
%  Example 3: Noisy objective function
%  Example 4: Noisy objective with user-provided noise estimates
%  Example 5: Periodic function
%  Example 6: Extended usage
%
%  Note: after installation, run 
%    bads('test') 
%  to check that everything is working correctly.
%
%  For any question, check out the FAQ: 
%  https://github.com/acerbilab/bads/wiki
%
%  See also BADS.

% Luigi Acerbi 2017-2022

disp('Running a number of examples usage for Bayesian Adaptive Direct Search (BADS).');
disp('Open ''bads_examples.m'' to see additional comments and instructions.');


%% Example 1: Basic usage

% Simple usage of BADS on Rosenbrock's banana function in 2D
% (see https://en.wikipedia.org/wiki/Rosenbrock_function).

% Screen display
fprintf('\n');
disp('*** Example 1: Basic usage');
disp('  Simple usage of BADS on <a href="https://en.wikipedia.org/wiki/Rosenbrock_function">Rosenbrock''s banana function</a> in 2D.');
disp('  Press any key to continue.'); fprintf('\n');
pause;

% We specify wide hard bounds and tighter plausible bounds that (hopefully) 
% contain the solution. Plausible bounds represent your best guess at 
% bounding the region where the solution might lie.

x0 = [0 0];                 % Starting point
lb = [-20 -20];             % Lower bounds
ub = [20 20];               % Upper bounds
plb = [-5 -5];              % Plausible lower bounds
pub = [5 5];                % Plausible upper bounds


% Run BADS, which returns the minimum X and its value FVAL.
[x,fval] = bads(@rosenbrocks,x0,lb,ub,plb,pub)

disp('The true global minimum is at X = [1,1], where FVAL = 0.');

% Note that BADS by default does not aim for extreme numerical precision 
% (e.g., beyond the 2nd or 3rd decimal place), since in realistic 
% model-fitting problems such a resolution is typically pointless.


%% Example 2: Non-bound constraints

fprintf('\n');
disp('*** Example 2: Non-bound constraints');
disp('  As before, but we force the input to stay in a circle with unit radius.');
disp('  BADS will complain because the plausible bounds are not specified explicitly.');
disp('  Press any key to continue.'); fprintf('\n');
pause;

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

% Run BADS with both bound and non-bound constraints
[x,fval] = bads(@rosenbrocks,x0,lb,ub,[],[],nonbcon)

% Alternatively, the following instructions would make BADS happier 
% (in a realistic model-fitting scenario, we recommend, whenever possible, 
% to specify plausible bounds which are tighter than the hard bounds).
%
% plb = lb; pub = ub;
% [x,fval] = bads(@rosenbrocks,x0,lb,ub,plb,pub,nonbcon)

disp('The true global minimum under these constraints is at X = [0.786,0.618], where FVAL = 0.046.');


%% Example 3: Noisy objective function

% We test BADS on a noisy function. For the purpose of this test, we
% manually add unit Gaussian noise to the 'sphere' (quadratic) function.

fprintf('\n');
disp('*** Example 3: Noisy objective function');
disp('  We test BADS on a noisy quadratic function with unit Gaussian noise.');
disp('  Press any key to continue.'); fprintf('\n');
pause;

% Define the noisy objective as a function handle
noisyfun = @(x) sum(x.^2,2) + randn(size(x,1),1);

% For a change, we start farther away from the solution
x0 = [-3 -3];

% This time to help the search we set tighter bounds
lb = [-5 -5];   ub = [5 5];
plb = [-2 -2];  pub = [2 2];

% We set several OPTIONS for the optimization (start with defaults)
options = bads('defaults');

% For this optimization, we explicitly tell BADS that the objective is
% noisy (it is not necessary, but it is a good habit); and also specify a 
% rough estimate for the value of the standard deviation of the noise in a 
% neighborhood of the solution.
options.UncertaintyHandling = true;
options.NoiseSize = 1;  % Optional, leave empty if unknown

% We also limit the number of function evaluations, knowing that this is a 
% simple example. Generally, BADS will tend to run for longer on noisy 
% problems to better explore the noisy landscape.
options.MaxFunEvals = 300;

% Finally, we tell BADS to re-evaluate the target at the returned solution 
% with 100 samples (10 by default, but since our function is inexpensive we 
% can use more evaluations). Note that this number counts towards the budget 
% of function evaluations.
options.NoiseFinalSamples = 100;

% Run BADS on the noisy function
[x,fval,exitflag,output] = bads(noisyfun,x0,lb,ub,plb,pub,[],options);
x
fval
disp(['The true, noiseless value of the function at X is ' num2str(sum(x.^2,2)) '.']);
disp('The true global minimum is at X = [0,0], where FVAL = 0.');

% FVAL in this case is an *estimate* of the function at X, obtained by
% averaging OPTIONS.NoiseFinalSamples function evaluations (default 10). 
% These values can be found in the OUTPUT structure (OUTPUT.yval), together 
% with additional information about the optimization.

disp('  Press any key to continue.'); fprintf('\n');
pause;

disp('The returned OUTPUT structure is:');
output

% Note that the fractional overhead of BADS reported in OUTPUT is astronomical.
% The reason is that the objective function we are using is analytical and 
% extremely fast, which is not what BADS is designed for. 
% In a realistic scenario, the objective function will be moderately costly
% (e.g., more than 0.1 s per function evaluation), and the fractional 
% overhead should be less than 1.

% For more information on optimizing noisy objective functions, see the
% BADS wiki: https://github.com/acerbilab/bads/wiki#noisy-objective-function


%% Example 4: Noisy objective with user-provided noise estimates

% Sometimes you may be able to estimate the noise associated with *each* 
% function evaluation, for example via bootstrap or other estimation methods.
% If you can do that, it is highly recommended you do so and tell BADS.

fprintf('\n');
disp('*** Example 4: Noisy function with user-provided noise estimates');
disp('  We test BADS on a noisy function for which we can estimate the noise.');
disp('  Press any key to continue.'); fprintf('\n');
pause;

% First, we inform BADS that we are specifying the target noise.
options = bads('defaults');         % Get default options
options.UncertaintyHandling = true; % Tell BADS that the objective is noisy
options.SpecifyTargetNoise = true;  % We are also specifying the noise
options.NoiseFinalSamples   = 100;  % More samples at the end 

% The noise standard deviation is returned as *second output* of the target
% function (check out the function "hetsphere.m").

x0 = [-3 -3];
lb = [-5 -5];   ub = [5 5];
plb = [-2 -2];  pub = [2 2];

% Run BADS on the noisy function
[x,fval] = bads(@hetsphere,x0,lb,ub,plb,pub,[],options)

disp(['The true, noiseless value of the function at X is ' num2str(sum(x.^2,2)) '.']);
disp('The true global minimum is at X = [0,0], where FVAL = 0.');
disp('Due to the elevated level of noise, we do not necessarily expect high precision in the solution.');


%% Example 5: Objective function with periodic dimensions

% We test BADS on a function with a subset of periodic dimensions.

fprintf('\n');
disp('*** Example 5: Objective function with periodic dimensions');
disp('  We test BADS on a function with some periodic inputs.');
disp('  Press any key to continue.'); fprintf('\n');
pause;

% This function is periodic along the third and fourth dimension, with
% periods respectively 4 and 2.
periodicfun = @(x) rosenbrocks(x(:,1:2)) + cos(x(:,3)*pi/2) + cos(x(:,4)*pi) + 2;

x0 = [-3 -3 -1 -1];
% We specify the periodic bounds via the hard bounds
lb = [-10 -5 -2 -1];
ub = [5 10 2 1];

plb = [-2 -2 -2 -1];
pub = [2 2 2 1];

options = bads('defaults');         % Get default options
options.PeriodicVars = [3 4];       % The 3rd and 4th variables are periodic

[x,fval] = bads(periodicfun,x0,lb,ub,plb,pub,[],options)
disp('The true global minimum is at X = [1,1,±2,±1], where FVAL = 0.');


%% Example 6: Extended usage

% Extended usage of BADS that shows some additional options.

fprintf('\n');
disp('*** Example 6: Extended usage');
disp('  Extended usage of BADS with additional options and no detailed display.');
disp('  Press any key to continue.'); fprintf('\n');
pause;

% Function handle for function with multiple input arguments (e.g., here
% we add a translation of the input; but more in general you could be 
% passing additional data to your objective function).
fun = @(x,mu) rosenbrocks(bsxfun(@plus, x, mu)); 

% This will translate the Rosenbrock fcn such that the global minimum is at zero
mu = ones(1,4);

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
options.Display             = 'final';  % Print only basic output ('off' turns off)
options.UncertaintyHandling = false;    % The objective is deterministic

% Custom output function, this one just prints the iteration number 
% (return FALSE to continue, TRUE to stop optimization)
options.OutputFcn = @(x,optimState,state) ...
    ~isfinite(fprintf('%s %d... ', state, optimState.iter));

% Run BADS, passing MU as additional (fixed) input argument for FUN
[x,fval,exitflag,output] = bads(fun,x0,lb,ub,plb,pub,[],options,mu);

% The following line of code would do the same using an anonymous function
% [x,fval,exitflag] = bads(@(x) fun(x,mu),x0,lb,ub,plb,pub,[],options);

x
fval
disp('The true global minimum is at X = [0,0,0,0], where FVAL = 0.');

exitflag
disp('EXITFLAG of 0 means that the maximum number of function evaluations has been reached.');
fprintf('\n');

disp('For this optimization we used the following OPTIONS:')
options
disp('Type ''help bads'' for additional documentation on BADS, or consult the <a href="https://github.com/acerbilab/bads">Github page</a> or <a href="https://github.com/acerbilab/bads/wiki">online FAQ</a>.');