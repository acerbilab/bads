function xs = searchSeries(x,gpstruct,LB,UB,optimState,options)
%SEARCHSERIES Search step by time series prediction (UNSUPPORTED).

if nargin < 1
    xs = 'srs';
    return;
end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

xsuccess = bsxfun(@rdivide,optimState.xsuccess,optimState.scale);
fsuccess = optimState.fsuccess;

D = size(x,2);
t = size(xsuccess,1);

% Random rotation matrix
[R,~] = qr(randn(D));
R = eye(D);

fmu = zeros(1,D);
fs = zeros(1,D);

xsuccess = xsuccess*R';

for i = 1:D
    index = max(1,t-10):t;
    
    stdy = sqrt(std(xsuccess(index,i)).^2 + optimState.TolMesh.^2);
    gptemp = gpinit(stdy,options,optimState); 
    gptemp.x = (index)';
    % gptemp.x = -fsuccess(index); 
    gptemp.y = xsuccess(index,i);
    gptemp.hyp(1).lik = log(MeshSize);
    gptemp.prior.lik{1}{2} = gptemp.hyp(1).lik; 
    
    if numel(index) < 3 || 1
        lastindex = 1:numel(index);
    else
        lastindex = ceil(numel(index)/2):numel(index);
    end
    if numel(gptemp.hyp(1).mean) == 1
        gptemp.hyp(1).mean = mean(gptemp.y(lastindex));   
        gptemp.prior.mean{1}{2} = gptemp.hyp(1).mean;
    else
        gptemp.hyp(1).mean = polyfit(gptemp.x(lastindex),gptemp.y(lastindex),1)';
        gptemp.prior.mean{1}{2} = gptemp.hyp(1).mean(1);        
        gptemp.prior.mean{2}{2} = gptemp.hyp(1).mean(2);        
    end
    
    bounds = unwrap2vec(gptemp.bounds);
    lb = bounds(1:2:end-1);
    ub = bounds(2:2:end);                
    
    prior = @(theta_) independent_prior(gptemp.prior,theta_);
    gptemp.inf = {@inference_with_prior, @exact_inference, prior};    

    % hyp0(1) = gptemp.hyp;
    %s = min(max(unwrap2vec(independent_prior(gptemp.prior)),lb),ub);   
    %hyp0(2) = rewrap(gptemp.hyp,s);
    
    % Quick-and-dirty grid optimization
    % ell = 0:0.1:2;
    ell = 1;
    sf = log(stdy) + (-2:0.5:2);
    theta = combvec(ell,sf)';
    nlZ = zeros(1,size(theta,1));
    for j = 1:size(theta,1)
        nlZ(j) = gp_optimizer(theta(j,:), gptemp.hyp, gptemp.inf, gptemp.mean, gptemp.cov, gptemp.lik, gptemp.x, gptemp.y);
    end
    [~,index] = min(nlZ);
    gptemp.hyp.cov = theta(index,:)';
    
    cov(i,:) = gptemp.hyp.cov(:)';
    
    % optoptions = optimset('TolFun',0.01,'TolX',1e-4,'MaxFunEval',100);
    % optoptions.Hessian = 'user-supplied';
    % gptemp.hyp = gpHyperOptimize(hyp0,gptemp,optoptions); 
    
%   [~,~,fmu(i),fs2] = gpgrad(0,gptemp);
    dt = exp(SearchFactor-1);
    [~,~,fmu(i),fs2] = gpgrad(t+dt,gptemp);
    fs(i) = sqrt(fs2);
end

% cov

fs = sqrt(fs.^2 + (MeshSize*options.PollMeshMultiplier^(-options.SearchGridNumber)).^2);


xs = bsxfun(@plus,bsxfun(@times,randn(options.Nsearch,D),fs),fmu);

xs = xs*R;

if 1
    figure(5);
    hold off;
    %plot(-fsuccess,xsuccess(:,1:end-1));
    plot(1:t,xsuccess*R);
    hold on;
    scatter((t+dt)*ones(size(fmu)),fmu);    
    %scatter(0*ones(size(fmu)),fmu);    
    drawnow;
    set(gcf,'Color','w');
    set(gca,'TickDir','out');
end

end

%--------------------------------------------------------------------------
function gptemp = gpinit(stdy,options,optimState)
%BGA_GPINIT

TolFun = options.TolFun;
TolMesh = optimState.TolMesh;

%% gp covariance function

% Base covariance function properties: smoothness and length scales
% gptemp.cov = {@isotropic_matern_covariance, 5};
gptemp.cov = {@isotropic_sqdexp_covariance};
gptemp.covbase = gptemp.cov;

ncov = str2double(feval(gptemp.cov{1}));
gptemp.hyp.cov = zeros(ncov,1);

% Covariance priors
for i = 1:ncov-1
     gptemp.prior.cov{i} = {@gaussian_prior, 1, 2^2};
     gptemp.bounds.cov{i} = [-5;Inf];
end
gptemp.hyp.cov = [log(3);log(stdy)];
gptemp.prior.cov{ncov} = {@gaussian_prior, log(stdy), 2^2};
gptemp.bounds.cov{ncov} = [log(TolFun);Inf];

if 0
    gptemp.cov = {@covSum,{gptemp.covbase,{@covPERiso,{@covMaterniso,Smoothness*2+1}}}};
    gptemp.hyp.cov = [gptemp.hyp.cov;zeros(ncov+1,1)];
    for i = ncov+(1:ncov)
         gptemp.prior.cov{i} = {@gaussian_prior,-1,2^2};
         gptemp.bounds.cov{i} = [-5;Inf];
    end
    gptemp.prior.cov{ncov*2+1} = {@gaussian_prior,0,2^2};    
    gptemp.bounds.cov{ncov*2+1} = [log(TolFun);Inf];
elseif 0
    gptemp.cov = {@covariance_product,{gptemp.covbase,@NNone_covariance}};
    gptemp.hyp.cov = [gptemp.hyp.cov;zeros(ncov,1)];
    for i = ncov+(1:ncov-1)
         gptemp.prior.cov{i} = {@gaussian_prior,-1,2^2};
         gptemp.bounds.cov{i} = [-5;Inf];
    end
    gptemp.prior.cov{ncov*2} = {@gaussian_prior,0,2^2};    
    gptemp.bounds.cov{ncov*2} = [log(TolFun);Inf];
end

%% gp likelihood function

gptemp.lik = @likGauss;   % Gaussian likelihood
likmu = log(TolMesh); liks2 = 1^2;
gptemp.hyp.lik = likmu;
gptemp.prior.lik = {{@delta_prior, likmu}};
% gptemp.prior.lik = {{@gaussian_prior, likmu, liks2}};

gptemp.hyp.lik = gptemp.hyp.lik(:);
gptemp.bounds.lik = {[log(TolMesh); 10]};

%% gp mean function
if 1
    gptemp.mean = @constant_mean;             % Constant mean function
    gptemp.hyp.mean = NaN;
    gptemp.prior.mean = {{@delta_prior, 0}};      % Fixed mean
else
    gptemp.mean = {@meanSum,{@meanLinear,@meanConst}};
    gptemp.hyp.mean = NaN(2,1);
    gptemp.prior.mean{1} = {@delta_prior, 0};      % Fixed prior
    gptemp.prior.mean{2} = {@delta_prior, 0};
end

%% gp sampling weight
gptemp.hypweight = 1;
gptemp.hypmean = [];

%% gp inference method
gpinf = @inference_with_prior;              % Exact inference
gptemp.infMethod = 'exact';

prior = @(theta_) independent_prior(gptemp.prior,theta_);
gptemp.inf = {@inference_with_prior, @exact_inference, prior};

% gptemp.inf = {@infPrior,gpinf,gptemp.prior};

gptemp = gpstruct_check(gptemp, 1);

end

%--------------------------------------------------------------------------
function nlZ = gp_optimizer(theta, hyp, inf, mean, cov, lik, x, y)
%GP_OPTIMIZER Wrapper function for GP optimization

hyp.cov = theta(:);
f = @() (feval(inf{:}, hyp, mean, cov, lik, x, y));
[~, nlZ] = f();

end