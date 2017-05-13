% WARPTEST

%% Example 1: Basic usage

% Simple usage of BADS on Rosenbrock's banana function in 2D
% (see https://en.wikipedia.org/wiki/Rosenbrock_function).
% 
% We specify wide hard bounds and tighter plausible bounds that (hopefully) 
% contain the solution. Plausible bounds represent your best guess at 
% bounding the region where the solution might lie.

nvars = 3;
x0 = 5*ones(1,nvars);                 % Starting point
lb = -10*ones(1,nvars);               % Lower bounds
ub = 6*ones(1,nvars);                % Upper bounds
plb = 2.9*ones(1,nvars);              % Plausible lower bounds
pub = 3*ones(1,nvars);               % Plausible upper bounds

options.InputWarping = 1;
options.NonlinearScaling = 0;

% Screen display
%fprintf('\n');
%display('*** Example 1: Basic usage');
%display('  Simple usage of BADS on <a href="https://en.wikipedia.org/wiki/Rosenbrock_function">Rosenbrock''s banana function</a> in 2D.');
%display('  Press any key to continue.'); fprintf('\n');
%pause;

fun = @rosenbrocks;
fun = @(x) sum(sqrt(abs(x)),2);
[x,fval] = bads(@(x) fun(exp(x)),x0,lb,ub,plb,pub,[],options)

display('The true global minimum is at X = [1,1], where FVAL = 0.');
