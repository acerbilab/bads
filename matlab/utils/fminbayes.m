function [x,fval,exitFlag,output] = fminbayes(fun,x0,lb,ub,depvars,options,varargin)

% TODO
% 1) Add dimension 3 (not sure, actually)
% 2) Add UncertaintyHandling and Method options
% 3) Add checks if UncertaintyHandling/Method are not specified:
%    - Evaluate f(x0) three times
%    - Evaluate SD of f and mean running time
%    - If the function is stochastic (SD > TolFun) turn UncertaintyHandling on
%    - If UncertaintyHandling is off
%       - If RunTime > 0.5 s use gp, otherwise splines
%       - If RunTime ~ 0.5 warn the user that she'd better specify the Method
% 4) Add the 'Greedy' option
%    - Stop as soon an improvement is detected
%    - Could add bias when sequentially building a grid
%       - Solution: randomly shuffle the order of evaluation
% 5) For gp with UncertaintyHandling:
%    - For each point, record mean and SD of the latent function
%    - Sample/refine new points using EI as acquisition function
%       - Compute EI on a (finer?) grid, then maximize it
%       - Or should we use LUCB? (see below)
%       - I think that EI here makes sense: be brave when exploring, be
%         conservative when deciding which point to keep
%    - Poll next point according to *Lowest Upper* Confidence Bound (LUCB)
%       - We want to be sure that the taken point is good!
%       - Play around with different values of kappa (e.g., 1, 2 or 3)
%    - Think whether you want to keep Elitism of 1 -- perhaps it is better
%      without (otherwise a weird gp estimate can screw up everything)
%    - At the moment we use a fixed grid before acquiring the next point 
%      but this might change to a proper Bayesian Optimization setting:
%       - Initialize the grid on a 2^D hypercube
%       - Sequentially acquire N(D) points (or until improvement if Greedy)
%          - Do it in a batch way (e.g, 10 at a time) if 1-by-1 is too slow
%       - If N becomes too large, merge the closest two points in a single
%         measurement:
%            - x* is is weighted average of x1,x2 (based on sigma_n) 
%            - y* is weighted average of y1,y2 (based on sigma_n)
%            - sigma_n(x*) is updated accordingly
%       - For now there is no information retained between different
%         iterations, think about that
%          - We could use the (say) 50 points nearest to the line/plane, 
%            with distances rescaled by gp length scale in that direction
%          - For gp length scale, use geometric mean of the samples
%            (it is biased downright, which works for us, better to
%             underestimate than overestimate the relevance of a point)
%          - We do not refit all hyper-parameters, only the relevant ones
%          - We could actually use ALL (or a lot of) the points by
%            projecting them back onto the grid with increased noise, and
%            then merging them until they reach a small number
% 6) Put prior on hyp.lik parameter to prevent numerical problems when
%    training gp
% 7) Default options.TolMesh is 1e-5. 
% 8) GPINIT should be an external function, can be user-specified.
% 9) Create DATA struct that contains all observations, their noise, etc.
% 10) Fix the way in which the mean of the gp is chosen; weigh each
%     observation according to its precision.

x0 = x0(:)';
nvars = length(x0);

if nargin < 3 || isempty(lb); lb = -Inf(size(x0)); end
if nargin < 4 || isempty(ub); ub = Inf(size(x0)); end
if nargin < 5; depvars = []; end
if nargin < 6; options = []; end

% N = [8 6];
fval = Inf;
yval = Inf;
yvalold = Inf;
fvalold = Inf;

defopts.Scale = [];
% defopts.xlabels = {'\sigma_a^{(1)}','\sigma_a^{(2)}','\sigma_b','\sigma_{prior}','p_{common}','\lambda'};
defopts.xlabels = [];
%defopts.xtrue = [log([1.9238 5.7849 3.8117 4.0861]) 0.4927 0.0105];
defopts.gpSamples = 20;
defopts.Smoothness = 3;
defopts.NoisyFun = 1;
defopts.NoiseKnown = false;
defopts.Display = 'iter';
defopts.N = [18 6];
defopts.UncertaintyHandling = false;

defopts.TolFun = 1e-4;
defopts.TolMesh = 1e-3;
defopts.TolAccept = 1;
defopts.TolSkip = 0.25;
defopts.TolGreedyAccept = Inf;
%defopts.MeshExpansion = exp(log(2)/nvars);
%defopts.MeshContraction = 1/exp(log(2)/nvars);
defopts.MeshExpansion = 2; %exp(nvars*1.5*log(2)/nvars);
defopts.MeshContraction = 0.5; % exp(-log(2)/nvars);
%defopts.MeshContraction = 0.5;

defoptnames = fieldnames(defopts);
for i = 1:length(defoptnames)
    defo = defoptnames{i};
    if ~isfield(options,defo) || isempty(options.(defo))
        options.(defo) = defopts.(defo);
    end
end

lb = lb(:)';
ub = ub(:)';

if isempty(options.Scale);
    options.Scale = logscale(lb,ub,mean(x0(:),2))';
end

if options.TolGreedyAccept < options.TolAccept
    warning('TOLGREEDYACCEPT in OPTIONS struct should be at least as large as TOLACCEPT. Modified accordingly.');
    options.TolGreedyAccept = options.TolAccept;
end 

InitialMeshSize = min(options.Scale, (ub-lb));
MeshSize = InitialMeshSize;
scale = 1;

if nargin < 6; funwrapper = fun;
else funwrapper = @(x_) fun(x_,varargin{:}); end

if isempty(depvars)
    dims = num2cell(1:nvars);
else    
    unidims = setdiff(1:nvars,depvars(:));
    dims = num2cell(unidims);
    for i = 1:size(depvars,1); dims{end+1} = depvars(i,:); end
end

if isempty(options.xlabels)
    for i = 1:nvars; options.xlabels{i} = ['x_' num2str(i)]; end
end

dimiter = 1;
iter = 1;

yval = funwrapper(x0);
funcCount = 1;

xmin = x0;
ymin = yval;

if strcmpi(options.Display,'iter')
    displayFormat = ' %5.0f     %8s       %5.0f    %12.6g    %12.6g   %12.6g    %16s\n';
    fprintf(' Iteration   Variables    f-count          f(x)          E[f(x)]       SD[f(x)]          MeshSize\n');
    fprintf(displayFormat, ...
        0, '0', funcCount, yval, fval, Inf, '');
end

% Pattern search vectors (orthonormal base)
searchvectors = eye(nvars);

% dims

while 1
    
    % x0
    d = dims{dimiter};
    
    % warning off;
    
    % Evaluate function at initial search pattern
    [x,w,vec] = initialsearchvectors(x0,lb,ub,MeshSize,searchvectors(d,:),options.N);
    tic;
    [y,yvar] = evaluatefunc(funwrapper,x,options.NoiseKnown && 0);    
    runtime = toc/size(x,1);    
    
    if options.NoisyFun
        % Build gp
        gpmean = mean(y);
        gpcov = log([std(w,[],1)';std(y)]);
        % if options.NoiseKnown; gplik = log(mean(yvar))/2; else gplik = []; end
        if options.NoiseKnown; gplik = log(options.TolFun); else gplik = []; end
        gpstruct = gpinit(w,gpmean,gpcov,gplik,options.Smoothness,options.TolFun);
        
        % Build grid
        if numel(d) == 1            
            vecgrid = searchvectors(d,:).*MeshSize;            
            wgrid = linspace(min(w), max(w),201)';
            % wgrid = [-logspace(log10(abs(min(w))),-4,100),0,logspace(-4,log10(abs(max(w))),100)]';
            xgrid = bsxfun(@plus, x0, bsxfun(@times, wgrid, vecgrid));
        else
            vec1 = searchvectors(d(1),:).*MeshSize;
            vec2 = searchvectors(d(2),:).*MeshSize;
            w1 = linspace(min(w(:,1)),max(w(:,1)),21)';
            w2 = linspace(min(w(:,2)),max(w(:,2)),21)';
            w1l = repmat(w1,[1,N(2)]); w1l = w1l(:);
            w2l = repmat(w2',[N(2),1]); w2l = w2l(:);
            wgrid = [w1l,w2l];
            xgrid = repmat(x0,[N(2)*N(2),1]);
            xgrid = xgrid + bsxfun(@times,w1l,vec1) + bsxfun(@times,w2l,vec2);
        end
        
        % Add x0 to grid if it was not included
        if isempty(findvector(x0,xgrid))
            xgrid = [xgrid;x0];
            wgrid = [wgrid;zeros(1,size(wgrid,2))];
        end
                
        for r = 1:options.N(numel(d))
              
            % Train GP on data
            gpstruct.x = w;
            gpstruct.y = y;
            gpstruct.xi = wgrid;
            
            % sort(gpstruct.x)'
            
            % Remove x from xi
            for i = 1:size(gpstruct.x,1)
                index = findvector(gpstruct.x(i,:),gpstruct.xi);
                if ~isempty(index); gpstruct.xi(index,:) = []; end
            end            
            
            %gpstruct.x 
            %gpstruct.y
                        
            gpstruct = gptrain(gpstruct,options.gpSamples,gpstruct.knownNoise);
            Nhyp = length(gpstruct.hyp);
            
            yval = min(y);
            tradeoff = 0;

            % GP prediction for the grid
            %tic
            mz = zeros(Nhyp,size(gpstruct.xi,1));
            y2s = zeros(Nhyp,size(gpstruct.xi,1));
            hypw = zeros(Nhyp,1);
            for i = 1:Nhyp
                [mz(i,:),~,~,y2s(i,:)] = gp(gpstruct.hyp(i),gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y,gpstruct.xi);
                %mz(i,:) = gp(gpstruct.hyp(i),gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y,gpstruct.xi);
                sz = sqrt(y2s);
                hypw(i) = gpstruct.hyp(i).w; 
            end
            gammaz = (yval + tradeoff - mz)./sz;
            z = -sz.*(gammaz.*(0.5*erfc(-gammaz/sqrt(2))) + exp(-0.5*(gammaz.^2))/sqrt(2*pi));
            mz = mean(mz,1);
            
            z = sum(bsxfun(@times,hypw,z),1);
            
            [bestEI,index] = min(z);
            wbest = gpstruct.xi(index,:);
            % mz = mz(index);
            %toc
            
            % Stop polling if best expected improvement is too low
            if abs(bestEI) < min(options.TolAccept,options.TolSkip*(scale^2))
                % bestEI
                break;
            end

            % if ~isempty(findvector(wbest,w)); break; end

            xsearch = x0 + sum(bsxfun(@times,wbest(:),vec),1);
            ysearch = evaluatefunc(funwrapper,xsearch);
            x = [x; xsearch];
            w = [w; wbest];
            y = [y; ysearch];

            % Minimize from best point
            %optoptions = optimset('Display','off');
            %[wnew,fval,exitflag,output] = fmincon(@gpfun,w0,[],[],[],[],w(1,:),w(end,:),[],optoptions);

            % output

            % Stop polling if improvement is larger than greedy threshold 
            if (ymin - min(y)) > options.TolGreedyAccept*(scale^2);
                [ymin, min(y),options.TolGreedyAccept*(scale^2)]
                break; 
            end
            
            % [yi,ys2i] = gpfunsd(wnew);
            
        end
        
        if ~options.UncertaintyHandling
            %xsearch = x0 + sum(bsxfun(@times,wnew(:),vec),1);
            %ysearch = evaluatefunc(funwrapper,xsearch);
            %xall = [x; xsearch];
            %yall = [y(:); ysearch];
            [yval,index] = min(y);
            wsearch = w(index,:);
            
            % index
            
            xsearch = x0 + sum(bsxfun(@times,wsearch(:),vec),1);
        end
        
        ysd = zeros(1,Nhyp);
        for i = 1:Nhyp
            ysd(i) = exp(gpstruct.hyp(i).lik*2);
        end
        ysd = sqrt(mean(ysd)).*ones(size(w,1),1);
        
    else
        
        [fval,index] = min(y);
        if numel(d) == 2
            y = reshape(y,[N2,N2])';
            w0 = [w1l(index),w2l(index)];
            minfunc = @interp2fun;
        else
            w0 = w(index);
            minfunc = @interp1fun;        
        end
        optoptions = optimset('Display','off');
        [wnew,fval,exitflag,output] = fmincon(minfunc,w0,[],[],[],[],w(1,:),w(end,:),[],optoptions);
        % wnew = w(index);
        mz = y;
        ysd = zeros([size(w,1),1]);        
    end
                
    
    funcCount = funcCount + numel(y);

    scaleold = scale;    
    if (ymin - yval) > options.TolAccept*(scale^2)
        scale = scale*options.MeshExpansion;
        if abs(wsearch) < 0.25; InitialMeshSize(d(1)) = InitialMeshSize(d(1))/options.MeshExpansion;
        elseif abs(wsearch) > 0.75; InitialMeshSize(d(1)) = InitialMeshSize(d(1))*options.MeshExpansion; end
    else
        scale = scale*options.MeshContraction;
    end

    % max(0,(ymin - yval))/scale^2
    
    if yval < ymin
        x0 = xsearch;
        ymin = yval;
    end
    
    [fi,ys2i] = gpfunsd(wsearch);
    
    if strcmpi(options.Display,'iter')
        fprintf(displayFormat, ...
            iter, numarray2str(dims{dimiter},'%g',',','',''), funcCount, ymin, fi, sqrt(ys2i), num2str(scaleold,'%.4g'));
%        fprintf(displayFormat, ...
%            iter, numarray2str(dims{dimiter},'%g',',','',''), funcCount, ymin, fval, sqrt(ys2i), numarray2str(MeshSize(dims{dimiter}),'%.4g',',','',''));
    end
    
    MeshSize = InitialMeshSize*scale;
        
    % Plot function
    [word,index] = sort(w);
    mz = gpfun(word);    
    
    plotFun(d,dimiter,dims,iter,word*scale,mz,y(index),ysd(index),wsearch*scale,lb,ub,options);
    
    if scale < options.TolMesh; break; end
    
    dimiter = max(1,mod(dimiter+1,length(dims)+1));
        
    if dimiter == 1
        iter = iter + 1;
    end
    
end

x = x0;
exitFlag = 0;
output = [];

    function yi = gpfun(wi)
        yi = zeros(length(gpstruct.hyp),size(wi,1));
        for ii = 1:length(gpstruct.hyp)
            [yi(ii,:)] = gp(gpstruct.hyp(ii),gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y,wi);
        end
        yi = mean(yi,1);
    end

    function [yi,ys2i] = gpfunsd(wi)
        yi = zeros(length(gpstruct.hyp),1);
        ys2i = zeros(length(gpstruct.hyp),1);
        for ii = 1:length(gpstruct.hyp)
            [yi(ii),~,~,ys2i(ii)] = gp(gpstruct.hyp(ii),gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y,wi);
        end
        yi = mean(yi,1);
        ys2i = mean(ys2i,1);
    end


    function yi = interp1fun(wi)
        yi = interp1(w,y,wi,'splines');
    end

    function yi = interp2fun(wi)
        yi = interp2(w1,w2,y,wi(1),wi(2),'splines');
    end

end

%--------------------------------------------------------------------------
function plotFun(D,dimiter,dims,iter,w,mz,y,ysd,wnew,lb,ub,options)
%PLOTFUN Plot figure

nrows = 2;
ncols = ceil(length(dims)/2);

switch numel(D)
    case 1
        subplot(nrows,ncols,dimiter);
        hold off;                
        plot(w,mz,'k','LineWidth',0.5); hold on;
        
        errorbar(w,y,ysd,'k','LineStyle','none');
        ylims = ylim(); % [min(mz),max(mz)]
        if ~isempty(options.xtrue)
            plot(options.xtrue(D)*[1 1],ylims,'r-','LineWidth',2);
        end
        plot(wnew*[1 1],ylims,'k-','LineWidth',1);
        xlim([lb(D),ub(D)]);
        xlabel(options.xlabels{D});
        
    case 2
        d = D(1);
        subplot(nrows,ncols,dimiter);
        hold off;
        %surf(w(:,1),w(:,2),y); hold on;
        %errorbar(w,y,ysd,'k','LineStyle','none');        
        zlims = zlim(); % [min(mz),max(mz)]
        if ~isempty(options.xtrue)
            plot3(options.xtrue(D(1)),options.xtrue(D(2)),1e5,'ro','LineWidth',2,'MarkerFaceColor','r');
        end
        colormap(1-bone);
        plot3(wnew(1),wnew(2),1e5,'ko','LineWidth',2,'MarkerFaceColor','k');
        view([0 90]);
        xlim([lb(D(1)),ub(D(1))]);
        ylim([lb(D(2)),ub(D(2))]);
        xlabel(options.xlabels{D(1)});
        ylabel(options.xlabels{D(2)});
end

set(gcf,'Color','w');
set(gca,'TickDir','out');
box off;
subplot(nrows,ncols,ceil(ncols/2)); title(['Iteration ' num2str(iter) ', D = ' num2str(dimiter)]);
for iRow = 1:nrows;
    subplot(nrows,ncols,1+(iRow-1)*ncols); ylabel('nLL');
end

drawnow;

end

%--------------------------------------------------------------------------
function gpstruct = gpinit(w,gpmean,gpcov,gplik,smoothness,tolfun)
%BGA_GPINIT Initialize gp struct for Bayesian genetic algorithm.

nvars = length(gpcov)-1;                    % Number of dimensions

%% gp mean function
gpstruct.meanfunc = @meanConst;             % Constant mean function
gpstruct.hyp.mean = gpmean;
gpstruct.prior.mean = {{@priorDelta}};      % Fixed mean

%% gp covariance function

switch smoothness
    case {0,1,2}                            % Matern-ARD covariance
        gpstruct.covfunc = {@covMaternard, smoothness*2+1};
    case 3
        gpstruct.covfunc = @covSEard;      % SE-ARD covariance
end
gpstruct.hyp.cov(1:nvars) = gpcov(1:nvars);
gpstruct.hyp.cov(nvars+1) = gpcov(nvars+1);
gpstruct.hyp.cov = gpstruct.hyp.cov(:);

% Weak smooth box prior on covariance ARD parameters
width = w(end,:)-w(1,:);
for i = 1:nvars
    gpstruct.prior.cov{i} = {@priorSmoothBox1,log(width(i))-3,log(width(i)),0.1};
%    gpstruct.prior.cov{i} = [];
end
gpstruct.prior.cov{nvars+1} = {@priorSmoothBox1,log(tolfun),log(200/width(i)),0.1};
%gpstruct.prior.cov{nvars+1} = []; % Flat prior on signal standard deviation

%% gp likelihood function

gpstruct.likfunc = @likGauss;   % Gaussian likelihood
if ~isempty(gplik)   % Known noise level
    gpstruct.prior.lik = {{@priorDelta}};
    gpstruct.hyp.lik = gplik;
    gpstruct.knownNoise = true;
else                % Unknown noise level
    likmu = log(1); liks2 = 1^2;
    gpstruct.hyp.lik = likmu;
    gpstruct.prior.lik = {{@priorGauss,likmu,liks2}};
    gpstruct.knownNoise = false;
end
gpstruct.hyp.lik = gpstruct.hyp.lik(:);

%% gp sampling weight
gpstruct.hyp.w = 1;        
gpstruct.hypmean = [];

gpstruct.lb = [log(min(width)/200),log(tolfun/max(width)),log(tolfun)];
gpstruct.ub = [log(max(width)),Inf,Inf];

%% gp inference method

if 1
    gpinf = @infExact;              % Exact inference
    gpstruct.infMethod = 'exact';
else
    gpstruct.infMethod = 'FITC_EP';
    n = 20;
    nu = fix(n/2); u = linspace(w(1,:),w(end,:),nu)';
    gpstruct.covfunc = {@covFITC, {gpstruct.covfunc}, u};
    gpinf = @infFITC_EP;          % also @infFITC_Laplace is possible    
end


gpstruct.inf = {@infPrior,gpinf,gpstruct.prior};

end

%--------------------------------------------------------------------------
function gpstruct = gptrain(gpstruct,Nsamples,NoiseKnown)
%GPTRAIN Train Gaussian Process hyper-parameters (optimize or sample).

switch lower(gpstruct.infMethod)
    case 'fitc_ep'
        samples = gpminimize(gpstruct.hyp,@gp,-300,0,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);
%  [a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, t, ones(n,1));
  
    case 'exact'
        ncalls = 0;
        if Nsamples == 0
            hyp = gpstruct.hyp;
            covlen = length(hyp.cov);
            if NoiseKnown; x0 = hyp.cov;
            else x0 = [hyp.cov;hyp.lik]; end
            %gpstruct.hyp = gpminimize(gpstruct.hyp,@gp,-100,0,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);
            samples = gpminimize(x0,@slice_wrapper,-300,0,gpstruct,-1);
        else
            covlen = length(gpstruct.hyp(1).cov);
            if length(gpstruct.hyp) == 1 || 1
                hyp = gpstruct.hyp(end);
                if NoiseKnown; x0 = hyp.cov;
                else x0 = [hyp.cov;hyp.lik]; end
                if ~isempty(gpstruct.hypmean); x0 = gpstruct.hypmean(:); end
                
                %tic
                %samples = gpminimize(x0,@slice_wrapper,-300,0,gpstruct,-1);
                %a = toc;                
                %tic                
                    uncopt = optimset('Display','off','GradObj','on','TolFun',1e-5,'TolX',1e-4);
                    %samples = fminunc(@(x_) slice_wrapper(x_,gpstruct,-1),x0,uncopt);
                    samples = fmincon(@(x_) slice_wrapper(x_,gpstruct,-1),x0,[],[],[],[],gpstruct.lb,gpstruct.ub,[],uncopt);
                %b = toc;
                %[a b]
                
                % hh = hessian(@(x_) slice_wrapper(x_,gpstruct,-1),samples);
                hyp.cov = samples(1:covlen,1);
                if ~NoiseKnown
                    hyp.lik = samples(covlen+1:end,1);
                    hyp.lik = max(hyp.lik,-5);
                end
                hyp.w = 1;
            else
                hyp = gpstruct.hyp(end);
            end
            
            % hyp = gpstruct.hyp;

            if NoiseKnown
                x0 = hyp.cov;
            else
                x0 = [hyp.cov;hyp.lik];
            end
            widths = 2*ones(size(x0));
                            
            if 1
            % vec = combvec([-1 0 1],[-1 0 1],[-1 0 1])';
                vec = combvec([-2 -1 0 1 2],[-2 -1 0 1 2],[-2 -1 0 1 2])';
                x = bsxfun(@plus,x0',vec);
            
                if 0
                    for d = 1:length(x0)
                        xx = linspace(x0(d)-widths(d),x0(d)+widths(d),21);
                        for i = 1:length(xx)
                            xt = x0; xt(d) = xx(i);
                            lZ(i) = slice_wrapper(xt,gpstruct);
                        end
                        subplot(1,length(x0),d);
                        plot(xx,lZ);
                        ylim([max(lZ)-10,max(lZ)]);
                        drawnow;
                    end
                end
                
                for i = 1:size(x,1)
                    lZ(i) = slice_wrapper(x(i,:),gpstruct);                    
                end
                hypw = exp(lZ)./sum(exp(lZ));
                gpstruct.hypmean = sum(bsxfun(@times,hypw(:),x),1);
                %hypw
                %x
                %gpstruct.hypmean
                                
                [hypword,index] = sort(hypw,'descend');
                samples = x(index(1:Nsamples),:)';
                hypword = hypword(1:Nsamples)/sum(hypword(1:Nsamples));
                
                cumw = 1 - cumsum(hypword);
                cumw(cumw < 0.01) = [];
                Nsamples = length(cumw);
                hypword(Nsamples+1:end) = [];
                hypword = hypword/sum(hypword);
                gpstruct.hyp(Nsamples+1:end) = [];
                
            else
                word = [];
                while 1
                    %try
                        ncalls = 0;
                        samples = slice_sample(Nsamples,Nsamples,2,@slice_wrapper,x0,widths,0,0,gpstruct);
                        %c = cov(samples');
                        %for i = 1:size(c,1)
                        %    for j = 1:size(c,2)
                        %        co(i,j) = c(i,j)./sqrt(c(i,i)*c(j,j));
                        %    end
                        %end
                        %co
                        break;
                    %catch
                    %    warning('Error in sampling... trying to recover.');
                    %    x0 = zeros(1,length(x0));
                    %end
                end
            end
            % samples
            % ncalls
        end
        
        for iSamp = 1:Nsamples
            gpstruct.hyp(iSamp) = hyp;
            gpstruct.hyp(iSamp).cov = samples(1:covlen,iSamp);
            if ~NoiseKnown
                gpstruct.hyp(iSamp).lik = samples(covlen+1:end,iSamp);
                gpstruct.hyp(iSamp).lik = max(gpstruct.hyp(iSamp).lik,-5);
            end
            if ~isempty(hypword); gpstruct.hyp(iSamp).w = hypword(iSamp); end
        end
                
end

    function [lZ,dlZ] = slice_wrapper(x,gpstruct,sig)
        ncalls = ncalls + 1;
        if nargin < 3; sig = 1; end
        hyp.cov = x(1:covlen);
        if ~NoiseKnown; hyp.lik = x(covlen+1:end); end
        hyp.lik = max(hyp.lik,-5);
        if nargout == 1
            nlZ = gp(hyp,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);            
        else
            [nlZ,ndlZstruct] = gp(hyp,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);
            dlZ = ndlZstruct.cov(:);
            if ~NoiseKnown; dlZ = [dlZ; ndlZstruct.lik(:)]; end
            dlZ = sig*(-dlZ);
        end
        lZ = sig*(-nlZ);
    end

    function [dlZ] = slice_wrapper_grad(x,gpstruct,sig)
        ncalls = ncalls + 1;
        if nargin < 3; sig = 1; end
        hyp.cov = x(1:covlen);
        if ~NoiseKnown; hyp.lik = x(covlen+1:end); end
        hyp.lik = max(hyp.lik,-5);
        [nlZ,ndlZstruct] = gp(hyp,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);
        dlZ = ndlZstruct.cov(:);
        if ~NoiseKnown; dlZ = [dlZ; ndlZstruct.lik(:)]; end
        dlZ = sig*(-dlZ);
    end


end

%--------------------------------------------------------------------------
function scale = logscale(LB,UB,meanX)
%LOGSCALE Used to determine the scaling factor for mesh

%   Copyright 2003-2007 The MathWorks, Inc.

%The following is a guideline for scaling variables.

%meanX, LB and UB could be (in)finite and (non)zero. Scaling will be different in each case. 
LowUpFiniteNonzero = (isfinite(LB) & abs(LB)>=eps) & (isfinite(UB)  & abs(UB)>=eps);
LowFiniteNonZero   = (isfinite(LB) & abs(LB)>=eps) & (~isfinite(UB) | abs(UB)<=eps);
UpFiniteNonZero    = (~isfinite(LB)| abs(LB)<=eps) & (isfinite(UB)  & abs(UB)>=eps);
InfiniteZero       = (~isfinite(LB)| abs(LB)<=eps) & (~isfinite(UB) | abs(UB)<=eps);
%Calculate LOG scale (Dennis & Schnabel)
lscale = zeros(length(meanX),1);
lscale(LowUpFiniteNonzero) = (log2(abs(LB(LowUpFiniteNonzero)))+log2(abs(UB(LowUpFiniteNonzero))))/2;
lscale(LowFiniteNonZero)   =  log2(abs(LB(LowFiniteNonZero)));
lscale(UpFiniteNonZero)    =  log2(abs(UB(UpFiniteNonZero)));

if abs(meanX(InfiniteZero))>=eps
lscale(InfiniteZero)       =  log2(abs(meanX(InfiniteZero)));
end
%Convert to normal scale.
scale = 2.^round(lscale);
end

%--------------------------------------------------------------------------
function [minw,maxw] = boundvec(v,x,lb,ub)
%BOUNDVEC Return bounds on coefficient for vector V starting from X
%  such that X + w*V is within the box defined by LB and UB
nvars = length(v);
maxw = Inf(1,nvars);
minw = -Inf(1,nvars);
pos = v > 0; neg = v < 0;
maxw(pos) = (ub(pos) - x(pos))./v(pos);
minw(neg) = (ub(neg) - x(neg))./v(neg);
minw(pos) = -(x(pos) - lb(pos))./v(pos);
maxw(neg) = -(x(neg) - lb(neg))./v(neg);
maxw = min(maxw); minw = max(minw);
end

%--------------------------------------------------------------------------
function [x,w,vec] = initialsearchvectors(x0,lb,ub,MeshSize,searchdir,N)

ndims = size(searchdir,1);

switch ndims
    case 1
        vec = searchdir.*MeshSize;
        [minw,maxw] = boundvec(vec,x0,lb,ub);

        if 0        
            w = linspace(max(minw,-1), min(maxw,1),N(1))';
            x = bsxfun(@plus, x0, bsxfun(@times, w, vec));
        else
            w = linspace(max(minw,-1), min(maxw,1),2)';
            x = bsxfun(@plus, x0, bsxfun(@times, w, vec));            
        end

    case 2
        vec1 = searchdir(1,:).*MeshSize;
        vec2 = searchdir(2,:).*MeshSize;
        [minw1,maxw1] = boundvec(vec1,x0,lb,ub);
        [minw2,maxw2] = boundvec(vec2,x0,lb,ub);
        w1 = linspace(max(minw1,-1),min(maxw1,1),N(2))';
        w2 = linspace(max(minw2,-1),min(maxw2,1),N(2))';
        w1l = repmat(w1,[1,N(2)]); w1l = w1l(:);
        w2l = repmat(w2',[N(2),1]); w2l = w2l(:);
        x = repmat(x0,[N(2)*N(2),1]);
        x = x + bsxfun(@times,w1l,vec1) + bsxfun(@times,w2l,vec2);
        w = [w1l,w2l];
        vec = [vec1;vec2];
        
    otherwise
        error('The size of D should be 1 or 2.');
end
end

%--------------------------------------------------------------------------
function [y,yvar] = evaluatefunc(fun,x,noise)
%EVALUATEFUNC Evaluate function at several points

if nargin < 3; noise = []; end
if isempty(noise) && nargout < 2; noise = 0; end

n = size(x,1);
y = NaN(n,1);
yvar = NaN(n,1);

if noise
    for iIter = 1:n
        [yt,yvart] = fun(x(iIter,:));
        y(iIter) = sum(yt); yvar(iIter) = sum(yvart);
    end
else
    for iIter = 1:n
        yt = fun(x(iIter,:));
        y(iIter) = sum(yt);
    end
end
end

%--------------------------------------------------------------------------
function index = findvector(x,y,tol)
%FINDVECTOR Returns position of vector X in array Y, within numerical
%   precision TOL (default 1e-12).

if nargin < 3 || isempty(tol); tol = 1e-12; end
if ~isvector(x); error('X should be a vector.'); end
if size(x,2) ~= size(y,2); error('X and Y should be arrays of the same length.'); end

index = find(mean(abs(bsxfun(@minus,y,x)),2) < tol);

end