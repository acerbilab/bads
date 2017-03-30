function [us,optimState] = searchMax(acq,u,gpstruct,LB,UB,optimState,options)

%SEARCHEI Search step by local maximization of expected improvement.

if nargin < 2 || isempty(u) || isempty(gpstruct)
    us = 'max';
    return;
end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

SearchMeshSize = options.PollMeshMultiplier.^optimState.SearchSizeInteger;
LB = LB(:) + SearchMeshSize/2;
UB = UB(:) - SearchMeshSize/2;
nvars = size(u,2);

%LB = max(LB, u(:) - 2*MeshSize*SearchFactor);
%UB = min(UB, u(:) + 2*MeshSize*SearchFactor);

[~,~,ftarget,fs2,~,post] = mygp(gpstruct.hyp(1),gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
    gpstruct.lik,gpstruct.x,gpstruct.y,u);
gpstruct.post = post;

ftarget = ftarget - optimState.sdlevel*sqrt(fs2 + options.TolFun^2);

if 0
    
    optoptions = optimset('Display','off','GradObj','on','DerivativeCheck','off',...
        'TolX',optimState.TolMesh,'TolFun',options.TolFun);
    
    idx = 1;
    funcCount = 0;
    while funcCount < options.Nsearch && idx <= 10
        x0 = x0 + 0.1*Scale*randn(size(x0));
        try
            optoptions.MaxFunEval = options.Nsearch - funcCount;
            % [xs(idx,:),~,~,output] = fmincon(@(x_) LowestUpperBound(x_,gpstruct,Scale),x0,[],[],[],[],LB,UB,[],optoptions);
            [xs(idx,:),~,~,output] = fmincon(@(x_) NegExpectedImprovement(x_,gpstruct),x0,[],[],[],[],LB,UB,[],optoptions);
            funcCount = funcCount + output.funcCount;
        catch
            warning('ah');
            xs(idx,:) = x0;
        end    
        idx = idx + 1;
    end

elseif 0
    
    prt = 1;

    smax = 10 + 5*nvars;            % Default maximum number of levels
    nf = options.Nsearch;           % Maximum number of fcn evaluations
    stop = 10 + 2*nvars;            % Stall iteration limit
    
    iinit.x0 = [0.5*(LB(:) + u(:)), 0.5*(LB(:) + u(:)) + 0.5*rand(nvars,1).*(UB(:)-LB(:)), 0.5*(UB(:) + u(:))];
    % iinit.x0 = [LB(:), x(:), UB(:)];
    iinit.l = 2*ones(nvars,1);
    iinit.L = 3*ones(nvars,1);

    data.acq = acq;
    data.target = ftarget;
    data.gpstruct = gpstruct;
    data.optimState = optimState;    
    
    [us,fbest,xmin,fmi,ncall,ncloc]=...
        mymcs('searchMaxAcq',data,LB(:),UB(:),prt,smax,nf);
%        mymcs('searchMaxAcq',data,LB(:),UB(:),prt,smax,nf,stop,iinit);
    
    us = us(:)';
    
else
        
    PLB = LB;
    PUB = UB;    
    data.acq = acq;
    data.target = ftarget;
    data.gpstruct = gpstruct;
    data.optimState = optimState;    
    cmaesopt.Display = 'iter';
    cmaesopt.TolFun = options.TolFun;
    cmaesopt.TolHistFun = options.TolFun/10;    
    [us,fbest] = cmaes_wrapper(@searchMaxAcq,u,LB,UB,PLB,PUB,cmaesopt,data);    
    
end