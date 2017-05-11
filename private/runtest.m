function failed = runtest()
%RUNTEST Test Bayesian Adaptive Direct Search (BADS).
%  RUNTEST executes a few runs of the BADS optimization algorithm to
%  check that it is installed correctly, and returns the number of failed
%  tests.
% 
%  See also BADS, BADS_EXAMPLES.

nvars = 3;                              % Number of dimensions
x0 = 4*ones(1,nvars);                  % Initial point
tolerr = [0.1 0.1 1];                   % Error tolerance

txt{1} = 'Test with deterministic function (ellipsoid)';
fprintf('%s.\n', txt{1});
fun = @(x) sum((x./(1:numel(x)).^2).^2);     % Objective function
[exitflag(1),err(1)] = testblock(fun,[],[],x0,0);

txt{2} = 'Test with deterministic function (sphere) and non-bound constraints';
fprintf('%s.\n', txt{2});
fun = @(x) sum(x.^2,2);
nonbcon = @(x) x(:,1) + x(:,2) < sqrt(2);     % Non-bound constraints
[exitflag(2),err(2)] = testblock(fun,[],nonbcon,x0,1);

txt{3} = 'Test with noisy function (noisy sphere)';
fprintf('%s.\n', txt{3});
fun = @(x) sum(x.^2) + randn();             % Noisy objective function
truefun = @(x) sum(x.^2);
[exitflag(3),err(3)] = testblock(fun,truefun,[],x0,0);

failed = 0;
fprintf('===========================================================================\n');
for i = 1:3
    fprintf('%s:', txt{i});
    if exitflag(i) >= 0 && err(i) < tolerr(i)
        fprintf('\tPASSED\n');
    else
        fprintf('\tFAILED\n');
        failed = failed + 1;
    end 
end
fprintf('===========================================================================\n');
fprintf('\n');

if failed == 0
    display('BADS is working correctly. See bads_examples.m for usage examples; check out the <a href="https://github.com/lacerbi/bads">BADS website</a>;'); 
    display('consult the <a href="https://github.com/lacerbi/bads/wiki">online FAQ</a>; or digit ''help bads'' for more information. Enjoy!');
else
    display('BADS is not working correctly. Please check the <a href="https://github.com/lacerbi/bads/wiki">online FAQ</a> for more information.');
end
    
    
% %% Noisy sphere (FMINSEARCH)
% 
% display('Comparison with FMINSEARCH (noisy sphere). Press any key to continue.')
% fun = @(x) sum(x.^2) + randn();             % Noisy objective function
% options = optimset('Display','iter');
% pause;
% [x,fval,exitflag,output] = fminsearch(fun,x0,options);
% display(['Final value (not-noisy): ' num2str(sum(x.^2),'%.3f') ' (true value: 0.0) with ' num2str(output.funcCount) ' fun evals.']);
% 
% %% Noisy sphere (GA)
% 
% if exist('ga.m','file')
%     display('Comparison with GA (noisy sphere). Press any key to continue.')
%     fun = @(x) sum(x.^2) + randn();             % Noisy objective function
%     options = gaoptimset('Display','iter','Generations',19);    
%     pause;
%     rng(0);
%     [x,fval,exitFlag,output] = ga(fun,nvars,[],[],[],[],LB,UB,[],[],options);
%     display(['Final value (not-noisy): ' num2str(sum(x.^2),'%.3f')  ' (true value: 0.0) with ' num2str(output.funccount) ' fun evals.']);
% end

end

%--------------------------------------------------------------------------
function [exitflag,err] = testblock(fun,truefun,nonbcon,x0,fmin)

nvars = numel(x0);
LB = -100*ones(1,nvars);                % Lower bound
UB = 100*ones(1,nvars);                 % Upper bound
PLB = -8*ones(1,nvars);                 % Plausible lower bound
PUB = 12*ones(1,nvars);                 % Plausible upper bound

options = bads('defaults');             % Default options
options.MaxFunEvals = 50;

[x,fval,exitflag,output] = bads(fun,x0,LB,UB,PLB,PUB,nonbcon,options);

if isempty(truefun)
    display(['Final value: ' num2str(fval,'%.3f') ' (true value: ' num2str(fmin) '), with ' num2str(output.funccount) ' fun evals.']);
    err = abs(fval - fmin);
else
    fval_true = truefun(x);
    display(['Final value (not-noisy): ' num2str(fval_true,'%.3f') ' (true value: ' num2str(fmin) ') with ' num2str(output.funccount) ' fun evals.']);
    err = abs(fval_true - fmin);
end
fprintf('\n');
end