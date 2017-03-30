function xs = searchOptim(x0,gpstruct,LB,UB,Scale,optimState,options)
%SEARCHOPTIM Search step by local maximization of expected improvement.

if nargin < 1
    xs = 'ei';
    return;
end

optoptions = optimset('Display','off','GradObj','on','DerivativeCheck','off',...
    'TolX',optimState.TolMesh,'TolFun',options.TolFun);

[ymui,ys2i,fmui,fs2i,~,post] = gp(gpstruct.hyp(1),gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
    gpstruct.lik,gpstruct.x,gpstruct.y,x0);
gpstruct.post = post;

idx = 1;
funcCount = 0;
while funcCount < options.Nsearch && idx <= 10
    x0 = x0 + 0.1*Scale*randn(size(x0));
    try
        optoptions.MaxFunEval = options.Nsearch - funcCount;
        % [xs(idx,:),~,~,output] = fmincon(@(x_) LowestUpperBound(x_,gpstruct,Scale),x0,[],[],[],[],LB,UB,[],optoptions);
        [xs(idx,:),~,~,output] = fmincon(@(x_) NegExpectedImprovement(x_,optimState.ftarget,gpstruct),x0,[],[],[],[],LB,UB,[],optoptions);
        funcCount = funcCount + output.funcCount;
    catch
        warning('ah');
        xs(idx,:) = x0;
    end    
    idx = idx + 1;
end

end

%--------------------------------------------------------------------------
function [y,dy] = NegExpectedImprovement(xi,target,gpstruct,scale)
%NEGEXPECTEDIMPROVEMENT Return NFEI

if ~isempty(gpstruct.x0)
    xi = bsxfun(@minus,xi,gpstruct.x0);
    gpstruct.x = bsxfun(@minus,gpstruct.x,gpstruct.x0);
end

Nhyp = length(gpstruct.hyp); 
hypw = zeros(Nhyp,1);

for i = 1:Nhyp; hypw(i) = gpstruct.hypweight(i); end
[ymu,ys2,fmu,fs2,dymu,dys2,dfmu,dfs2] = gpgrad(xi,gpstruct,'central');
fs = sqrt(fs2);
ys = sqrt(ys2);

gammaz = (target - fmu)./fs;
fpi = 0.5*erfc(-gammaz/sqrt(2));    % Probability of improvement
y = -fs.*(gammaz.*fpi + exp(-0.5*gammaz.^2)/sqrt(2*pi));
y = sum(bsxfun(@times,hypw,y),1);

if any(isnan(y) | isinf(y))
    y = 0;
    dy = zeros(1,size(xi,2));
    return;
end

dfs = 0.5*dfs2./fs;
dgammaz = -(dfmu.*fs + (target - fmu).*dfs)./fs2;
dfpi = -0.5*dgammaz/sqrt(2)*(-2*exp(-gammaz.^2/2)/sqrt(pi));

dy = -(dfs.*gammaz.*fpi + dgammaz.*fs.*fpi + dfpi.*gammaz.*fs) ...
    + (fs.*gammaz.*dgammaz - dfs).*exp(-0.5*gammaz.^2)/sqrt(2*pi);
dy = sum(bsxfun(@times,hypw,dy),1);

end

%--------------------------------------------------------------------------
function [y,dy] = LowestUpperBound(xi,gpstruct,scale,kappa)
%LOWESTUPPERBOUND Return LUB

if nargin < 4 || isempty(kappa); kappa = 1; end

if ~isempty(gpstruct.x0)
    xi = bsxfun(@minus,xi,gpstruct.x0);
    gpstruct.x = bsxfun(@minus,gpstruct.x,gpstruct.x0);
end

Nhyp = length(gpstruct.hyp); 
hypw = zeros(Nhyp,1);

for i = 1:Nhyp; hypw(i) = gpstruct.hypweight(i); end
if nargout > 1
    [ymu,ys2,fmu,fs2,dymu,dys2,dfmu,dfs2] = gpgrad(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2] = dgp(xi,gpstruct);    
end
fs = sqrt(fs2);
% ys = sqrt(ys2);

y = fmu + kappa.*fs;
y = sum(bsxfun(@times,hypw,y),1);

if any(isnan(y) | isinf(y))
    y = 0;
    dy = zeros(1,size(xi,2));
    return;
end

if nargout > 1
    dy = dfmu + 0.5*kappa.*dfs2./fs;
    dy = sum(bsxfun(@times,hypw,dy),1);
    if any(isnan(dy) | isinf(dy)); dy = zeros(1,size(xi,2)); end    
end

end