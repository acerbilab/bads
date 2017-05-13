function inpwarp = warpinput(optimState,options)

nvars = size(optimState.U,2);

% optimState.Y = sqrt(optimState.Y - min(optimState.Y));

% Define GP
gplik = [];
covtype = 'matern5';
covextras = 1;
gpstruct = gpdefBads(covtype,covextras,nvars,gplik,optimState,options);
u0 = optimState.U(1,:);
[gpstruct,exitflag] = gpTrainingSet(gpstruct,'all',u0,[],optimState,options,1);

% Initialize input warping
gpstruct.inpwarp.type = 1;
gpstruct.inpwarp.params = zeros(1,nvars*2);
gpstruct.inpwarp.prior.mu = zeros(1,nvars*2);
gpstruct.inpwarp.prior.sigma = 0.5*ones(1,nvars*2);
gpstruct.inpwarp.bounds.lb = -3*ones(1,nvars*2);
gpstruct.inpwarp.bounds.ub = 3*ones(1,nvars*2);

% Store lower/upper bounds
gpstruct.LB = optimState.LB;
gpstruct.UB = optimState.UB;

optoptions.Display = 'off';
optoptions.MaxFunEvals = min(50*nvars,300);

xbase = gpstruct.x;

for i = 1:3
    theta0 = gpstruct.inpwarp.params;
    gpstruct.x = xbase;
    try
        theta = fmincon(@(x) warp_optimizer(x, gpstruct),theta0,[],[],[],[],gpstruct.inpwarp.bounds.lb,gpstruct.inpwarp.bounds.ub,[],optoptions);
        gpstruct.inpwarp.params = theta;
    catch
        % Optimization failed, keep old warping
    end
    
    [gpstruct,exitflag] = gpTrainingSet(gpstruct,'all',u0,[],optimState,options,1);    
end

inpwarp = gpstruct.inpwarp;
gpstruct.inpwarp.params

end

%--------------------------------------------------------------------------
function nlZ = warp_optimizer(theta, gpstruct)
%WARP_OPTIMIZER Wrapper function for optimization of GP input warping

LB = gpstruct.LB;
UB = gpstruct.UB;

% Kumaraswmi parameters
a = exp(theta(1:2:end));
b = exp(theta(2:2:end));

% Warp inputs
xprime = uWarp(gpstruct.x,a,b,LB,UB);

% Compute marginal likelihood
try
    [~, nlZ] = feval(gpstruct.inf{:}, gpstruct.hyp, gpstruct.mean, gpstruct.cov, gpstruct.lik, xprime, gpstruct.y);
    
    % Add (negative) prior to negative marginal likelihood
    mu = gpstruct.inpwarp.prior.mu;
    sigma = gpstruct.inpwarp.prior.sigma;
    nlZ = nlZ + 0.5*sum(((theta - mu)./sigma).^2) + 0.5*sum(log(2*pi*sigma.^2));    
catch
    % Optimization failed
    nlZ = NaN;
end


end