function landscapeplot(fun,x0,LB,UB,scale,gpstruct,vec,n)

D = numel(x0);

if nargin < 5 || isempty(scale); scale = Inf; end
if nargin < 6; gpstruct = []; end
if nargin < 7 || isempty(vec); vec = eye(D); end
if nargin < 8 || isempty(n); n = 13; end

fval = zeros(1,n);
%lb = max(lb, x0-scale);
%ub = min(ub, x0+scale);
nd = size(vec,1);   % Number of explored landscape dimensions

% Random rotation matrix
%[R,~] = qr(randn(D));
%R = eye(D);
%vec = R*eye(D);


for i = 1:nd
    x{i} = repmat(x0,[n 1]);
    h{i} = scale*linspace(-1,1,n)';
    x{i} = bsxfun(@plus, x{i}, bsxfun(@times, h{i}, vec(i,:)));
    
    % i-th axis label
    if sum(vec(i,:) ~= 0) == 1
        % Axis vector
        xstr{i} = ['x_{' num2str(find(vec(i,:),1)) '}'];
    else
        % Generic vector
        xstr{i} = ['v_{' num2str(i) '}'];
    end
end

if ~isempty(gpstruct)
    gpxi = [];
    for i = 1:numel(x); gpxi = [gpxi; x{i}]; end
        
    %[ymu,ys2,fmu,fs2] = gp(gpstruct.hyp(1),gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
    %    gpstruct.lik,gpstruct.x,gpstruct.y,gpxi);    
    % fs = sqrt(fs2);
    
    [nfei,fmu,fs,ym,ys,fpi] = NegExpectedImprovement(gpxi,gpstruct);
    
end

for i = 1:nd
    if nd >= 4
        subplot(2,ceil(nd/2),i);
    elseif nd > 1
        subplot(1,nd,i);
    end
    
    for j = 1:n
        fval(j) = fun(x{i}(j,:));
    end
    miny = min(fval);
    fmintxt = ['f_{min} = ' num2str(miny,'%.4g')];

        
    hold off;
    if ~isempty(gpstruct)
        idx = (i-1)*n + (1:n);
        gpy = fmu(idx);
        gpsd = fs(idx);        
        ei = -nfei(idx);
        errorbar(h{i}, gpy, gpsd);
        hold on;        
        miny = min([miny;fmu(:)-fs(:)]);
        
        % GP Mean
        [~,idx] = min(gpy);
        fmintxt = [fmintxt ', E[f] = ' num2str(fval(idx),'%.4g')];
        
        % Expected improvement (EI)
        ei = ei*10;
        plot(h{i},ei+miny,'k','LineWidth',1);                
        [~,idx] = max(ei);
        fmintxt = [fmintxt ', f_{EI} = ' num2str(fval(idx),'%.4g')];
        
        % Lower upper bound (LUB)
        lub = gpy + gpsd;
        [~,idx] = min(lub);
        fmintxt = [fmintxt ', f_{LUB} = ' num2str(fval(idx),'%.4g')];
        
        if isfield(gpstruct,'nonlinf')
            fval(fval > gpstruct.nonlinmu) = ...
                gpstruct.nonlinf(fval(fval > gpstruct.nonlinmu),gpstruct.nonlinmu,gpstruct.deltay);
        end
    else
        fmu = []; fs = [];
    end
    
    plot(h{i}, fval, 'k-', 'LineWidth', 2); hold on;
    % plot([0 0],[miny,max([fmu(:) + fs(:);fval(:)])],'k--','LineWidth',0.5);
        
    box off;
    set(gca,'TickDir','out');
    xlabel(xstr{i});
    ylabel('f(x)');
    
    title(fmintxt);
    %ht = text(0, 0, ['f_{min} = ' num2str(min(fval),'%.4g')]);
    %set(ht, 'Units', 'normalized');
end

set(gcf,'Color','w');

end

%--------------------------------------------------------------------------
function [nfei,fm,fs,ym,ys,fpi] = NegExpectedImprovement(xi,gpstruct,eta)
%NEGEXPECTEDIMPROVEMENT Return NFEI

if nargin < 3 || isempty(eta); eta = 0; end

[ym,ys2,fm,fs2,hypw] = gppred(xi,gpstruct);
ys = sqrt(ys2);
fs = sqrt(fs2);

gammaz = (gpstruct.ftarget + eta - fm)./fs;
fpi = 0.5*erfc(-gammaz/sqrt(2));    % Probability of improvement
nfei = -fs.*(gammaz.*fpi + exp(-0.5*(gammaz.^2))/sqrt(2*pi));
nfei = sum(bsxfun(@times,hypw,nfei),1);
% ym = fm;

end

%--------------------------------------------------------------------------
function [lub,fm,fs,ym,ys,fpi] = LowestUpperBound(xi,gpstruct,kappa)
%LOWESTUPPERBOUND Return LUB

if nargin < 3 || isempty(kappa); kappa = 1; end

if isempty(gpstruct.x0)
    xx = gpstruct.x;
else
    xi = bsxfun(@minus,xi,gpstruct.x0);
    xx = bsxfun(@minus,gpstruct.x,gpstruct.x0);
end

Nhyp = length(gpstruct.hyp);
 
fm = zeros(Nhyp,size(xi,1));        % Function mean
fs = zeros(Nhyp,size(xi,1));        % Function SD
ys = zeros(Nhyp,size(xi,1));        % Observed value SD
hypw = zeros(Nhyp,1);

for i = 1:Nhyp
    try
        [fm(i,:),y2s,~,f2s] = gp(gpstruct.hyp(i),gpstruct.inf,gpstruct.mean,gpstruct.cov,gpstruct.lik,xx,gpstruct.y,xi);
    catch
        fm(i,:) = NaN*size(xi,1);
        y2s = NaN*size(xi,1);
        f2s = NaN*size(xi,1);
    end
    fs(i,:) = sqrt(f2s);
    ys(i,:) = sqrt(y2s);
    hypw(i) = gpstruct.hypweight(i);
end

ym = fm;
lub = fm + kappa.*fs;
lub = sum(bsxfun(@times,hypw,lub),1);  


end
