function [ymu,ys2,fmu,fs2,dymu,dys2,dfmu,dfs2] = gpgrad(x0,gpstruct,method,dx)
%GPGRAD Gradient of Gaussian process mean prediction

if nargin < 3 || isempty(method); method = 'central'; end
if nargin < 4 || isempty(dx); dx = 1e-6; end

D = size(x0,2);
Nhyp = length(gpstruct.hyp);

% Initialize variables
ymu = zeros(Nhyp,1);
fmu = zeros(Nhyp,1);
ys2 = zeros(Nhyp,1);
fs2 = zeros(Nhyp,1);

if nargout > 4
    dymu = zeros(Nhyp,D);
    dfmu = zeros(Nhyp,D);
    dys2 = zeros(Nhyp,D);
    dfs2 = zeros(Nhyp,D);
end

% Compute predictions and derivatives for each hyperparameter sample
for i = 1:Nhyp
    
    try
        if nargout > 4
            [dy,y0] = fgrad(@(xi_) gpfunc(xi_,gpstruct,gpstruct.hyp(i)),x0,method,'Vectorized','Step',dx); 
        else
            y0 = gpfunc(x0,gpstruct,gpstruct.hyp(i));
        end
    catch
        y0 = NaN(1,4);
        if nargout > 4; dy = NaN(4,D); end
        warning('bads:gpGradFail', 'Error in computing the GP derivative.');
    end
    
    ymu(i) = y0(1); 
    ys2(i) = y0(2);
    fmu(i) = y0(3);
    fs2(i) = y0(4);
    
    if nargout > 4
        dymu(i,:) = dy(1,:);
        dfmu(i,:) = dy(2,:);
        dys2(i,:) = dy(3,:);
        dfs2(i,:) = dy(4,:);
    end
        
        
end

%--------------------------------------------------------------------------

function yy = gpfunc(xi,gpstruct,hyp)
    if gpstruct.gpmlext
        % Pre-computed Laplace approximation at MAP
        if isfield(gpstruct,'hypSigma') && ~isempty(gpstruct.hypSigma)
            hyp.Sigma = gpstruct.hypSigma;
        end
        
        if isfield(gpstruct,'post') && ~isempty(gpstruct.post)
        %    [ymuib,ys2ib,fmuib,fs2ib] = mygp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        %gpstruct.lik,gpstruct.x,gpstruct.post,xi);
    
            [ymui,ys2i,fmui,fs2i] = mgp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,gpstruct.x,gpstruct.post,xi);
    
            %[fmuib(:) - fmui(:), sqrt(fs2ib(:)) - sqrt(fs2i(:))]
    
        else
            [ymui,ys2i,fmui,fs2i] = mgp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,gpstruct.x,gpstruct.y,xi);
        end
        yy = [ymui,ys2i,fmui,fs2i];
    else
        
        if isfield(gpstruct,'post') && ~isempty(gpstruct.post)
            [ymui,ys2i,fmui,fs2i] = mygp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,gpstruct.x,gpstruct.post,xi);
        else
            [ymui,ys2i,fmui,fs2i] = mygp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,gpstruct.x,gpstruct.y,xi);
        end
        
        if 1
            % [~, ymui, ys2i] = feval(gpstruct.lik{:}, hyp.lik, [], fmui, fs2i);
            
            lik = hyp.lik; lik(end) = -Inf;
            % fmuiold = fmui; fs2iold = fs2i;
            [~, fmui, fs2i] = feval(gpstruct.lik{:}, lik, [], fmui, fs2i);
        end
        
        yy = [ymui,ys2i,fmui,fs2i];
        % [ymui,sqrt(ys2i),fmui,sqrt(fs2i),fmuiold,sqrt(fs2iold)]
    end
    
    