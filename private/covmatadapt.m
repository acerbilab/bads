function optimState = covmatadapt(u,LB,UB,gpstruct,optimState,options)
%COVMATADAPT Covariance matrix adaptation (CMA).

% Estimation of local covariance structure (currently unused).

nvars = size(u,2);

cmasearch = searchCMA(1, ...
    u, ...
    gpstruct, ...
    LB, ...
    UB, ...
    optimState, ...
    options);
%toc

% optimState.C
% cmasearch

% Enforce periodicity
cmasearch = periodCheck(cmasearch,LB,UB,optimState);

% Force candidate points on search grid
cmasearch = force2grid(cmasearch, optimState);

% Remove already evaluated or unfeasible points from search set
cmasearch = uCheck(cmasearch,optimState.TolMesh,optimState,1);                                    

if size(cmasearch,1) > 2*nvars                     
    % Evaluate acquisition function on CMA search set
    cmaz = feval(options.SearchAcqFcn{:},cmasearch,optimState.ftarget,gpstruct,optimState,0);                
    cmaidx = isfinite(cmaz);
    cmaz = cmaz(cmaidx);
    if numel(cmaz) > 2*nvars
        cmasearch = cmasearch(cmaidx,:);

        % Compute vector weights
        mu = 0.5*size(cmasearch,1);
        weights = zeros(1,1,floor(mu));
        weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
        weights = weights./sum(weights);
        mueff = 1/sum(weights.^2);

        % Compute top vectors
        [~,cmaindex] = sort(cmaz,'ascend');

        % Compute weighted covariance matrix wrt X0
        Ubest = cmasearch(cmaindex(1:floor(mu)),:);
        Binv = ucov(Ubest,u,weights,optimState)./(optimState.meshsize*optimState.searchfactor).^2;
        
        Uworst = cmasearch(cmaindex(end:-1:end-floor(mu)+1),:);
        Uworst = [];
        
        if isempty(Uworst)
            amu = 1;
            cmu = min(0.2, amu*(mueff-2+1/mueff)/((nvars+2)^2 + amu*mueff/2));
            Binv = (1-cmu)*optimState.Binv + cmu*Binv;
        else
            Binv = Binv - ucov(Uworst,u,weights,optimState)./(optimState.meshsize*optimState.searchfactor).^2;
            cmu = 2/(nvars + sqrt(2))^2;
            beta = (4*mueff - 2)/((nvars+12)^2 + 4*mueff);
            Binv = (1-cmu)*optimState.Binv + beta*Binv;            
        end
        Binv = triu(Binv) + triu(Binv,1)'; % enforce symmetry
        try
            [V,D] = eig(optimState.Binv);           % eigen decomposition, B==normalized eigenvectors
            lambda = diag(D);                        
            %if isempty(lambda) || (min(lambda) < 0 && abs(min(lambda)/max(lambda)) > 1e-12)
            %    optimState.Binv = eye(nvars);
            %    optimState.C = eye(nvars);
            %    display('reset')

            lambdared = max(lambda,max(lambda)*1e-12);
            lambdared = sqrt(lambdared/sum(lambdared));
            lambdared = min(max(lambdared, optimState.searchmeshsize), max((optimState.UB-optimState.LB)./optimState.scale));
            lambdared = lambdared/sqrt(sum(lambdared.^2));
            optimState.Binv = Binv;
            optimState.C = real(V)*diag(lambdared);            
            % optimState.C
        catch
            % Failed, do nothing

        end
    end
end
