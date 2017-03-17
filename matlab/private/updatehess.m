function updatehess(uold,fvalold,fsdold,unew,fvalnew,fsdnew,optimState,gpstruct,options)
%UPDATEHESS Update estimated curvature (Hessian) at incumbent. (UNSUPPORTED)

nvars = size(unew,2);
    
% Compute derivatives at old and new point
g1 = fgrad(@(x_) gppred(x_,gpstruct), unew,'central','step',optimState.searchmeshsize);
g0 = fgrad(@(x_) gppred(x_,gpstruct), uold,'central','step',optimState.searchmeshsize);
g1 = g1./optimState.scale; g0 = g0./optimState.scale;
optimState.grad = g1(:); % Gradient at new incumbent
y = g1 - g0;
s = (unew - uold)./optimState.scale;
% sqrt(sum(s.*s))

if any(optimState.periodicvars)
    % Compute vector difference as the closest distance along periodic dimensions    
    % dual = optimState.UB - optimState.LB - abs(s);   % Opposite distance
    dual = (optimState.UB - optimState.LB)./optimState.scale - abs(s);   % Opposite distance
    for dim = find(optimState.periodicvars)
        if dual(dim) < abs(s(dim))      % Opposite distance is shorter
            s(dim) = -sign(s(dim))*dual(dim);
        end
    end
end

Binv = optimState.Binv;
switch options.HessianMethod
    case {'bfsg','bfgs','bgfs','bgsf','bsfg','bsgf'} % I never get it right
        if y*s' > 0
            Binv = Binv + (s*y' + y*Binv*y')*(s'*s)/(s*y')^2 - (Binv*y'*s + s'*y*Binv)/(s*y');
        else
            optimState.violations = optimState.violations + 1;
            % if optimState.violations < 1; return; end
            display(['negative ' num2str(optimState.violations)]);
            Binv = diag(gpstruct.lenscale.^2);
            optimState.violations = 0;
        end
        % Binv = 0.95*Binv + 0.05*eye(nvars) + (s*y' + y*Binv*y')*(s'*s)/(s*y')^2 - (Binv*y'*s + s'*y*Binv)/(s*y');
    case 'naive'
        c1 = 2/nvars^2; c2 = 0; c0 = c1+c2; c0 = 0;
        % gn = g1/sqrt(g1*g1');
        % Binv = Binv*(1-c0) + c1*(s'*s)/optimState.meshsize^2;
        Binv = Binv*(1-c0) + c1*(s'*s)/(s*s');
        if c2 > 0
            % Binv = Binv + c2*(gn'*gn)*((s*s')/optimState.meshsize)^2;            
            Binv = Binv + c2*diag(gpstruct.pollscale.^2); % /sqrt(mean(gpstruct.pollscale.^2));            
        end
        % s/optimState.meshsize
        % gn*optimState.meshsize
    case 'hessian'
        B = fhess(@(xi_) gppred(xi_,gpstruct), unew, [], 'step', optimState.searchmeshsize);
        c1 = 0;
        Binv = (1-c1)*Binv + c1*inv(B);
    case 'neighborhood'

        if 1
            optionstemp = options;
            optionstemp.Ndata = 2^10;
            optionstemp.MinNdata = 4 + floor(3*log(nvars));
            optionstemp.BufferNdata = Inf;
            optionstemp.gpRadius = 2;
            gptemp = gpTrainingSet(gpstruct, ...
                'neighborhood', ...
                uold, ...
                [], ...
                optimState, ...
                optionstemp, ...
                0);
            muratio = 0.25;
        else
            gptemp = gpstruct;
            muratio = 0.5;
        end
        X = gptemp.x;
        Y = gptemp.y;

        % Small jitter added to each direction
        jit = optimState.searchmeshsize;

        % Compute vector weights
        mu = floor(muratio*size(X,1));
        mu = max([mu, (4 + floor(3*log(nvars)))/2, floor(size(X,1)/2)]);
        weights = zeros(1,1,floor(mu));
        weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
        weights = weights./sum(weights);

        % Compute top vectors
        [~,index] = sort(Y,'ascend');

        % Compute weighted covariance matrix wrt X0
        u = uold;
        Ubest = X(index(1:floor(mu)),:);
        C = ucov(Ubest,u,weights,optimState)./optimState.meshsize^2;

        mueff = 1/sum(weights.^2);
        amu = 2; c1 = 0;
        % c1 = 2/((nvars+1.3)^2 + mueff); % Doesn't seem to improve --
        % might try to implement the whole path thing
        cmu = min(1-c1, amu*(mueff-2+1/mueff)/((nvars+2)^2 + amu*mueff/2));

        Binv = (1-cmu-c1)*Binv + c1*(s'*s)/optimState.meshsize^2 + cmu*C;

    otherwise
        error('Unknown Hessian update method.');
end

optimState.Binv = Binv;

if mod(optimState.iter,2) == 0 && options.HessianAlternate
    optimState.C = diag(gpstruct.pollscale./sqrt(sum(gpstruct.pollscale.^2)));
else        
    try
        [V,D] = eig(Binv);
        lambda = real(diag(D));
        % lambda(:)'
    catch
        lambda = [];
    end
    if isempty(lambda) || (min(lambda) < 0 && abs(min(lambda)/max(lambda)) > 1e-14)
        optimState.Binv = eye(nvars);
        optimState.B = eye(nvars);
        optimState.C = eye(nvars);
        display('reset')
    else
        lambda = max(lambda, max(lambda)*1e-14);
        optimState.B = real(V)*diag(1./lambda)*real(V)';
        lambdared = sqrt(lambda/sum(lambda));   
        lambdared = min(max(lambdared, optimState.searchmeshsize), max((optimState.UB-optimState.LB)./optimState.scale));
        % optimState.B = real(V)*diag(1./lambdared.^2)*real(V)';
        lambdared = lambdared/sqrt(sum(lambdared.^2));
        optimState.C = real(V)*diag(lambdared);
        optimState.Cres = optimState.C/sqrt(sum(1./lambdared));

        if any(~isreal(optimState.C(:)))
            optimState.Binv = eye(nvars);
            optimState.C = eye(nvars);
            display('unreal')
        end

        % log10(lambda)
    end
end