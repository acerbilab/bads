function [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,method,dx)
%GPPRED GP prediction

if nargin < 3; method = []; end
if nargin < 4; dx = []; end

n = size(xi,1);                 % Number of test points
Nhyp = numel(gpstruct.hyp);     % Number of hyper-parameter samples

% Requesting gp gradient?
compute_grad = nargout > 5;

% Rotate test input?
rotategp_flag = isfield(gpstruct,'C') && ~isempty(gpstruct.C); 
if rotategp_flag
    xi = xi*gpstruct.C';
end

if compute_grad
    [ymu,ys2,fmu,fs2,dymu,dys2,dfmu,dfs2] = gpgrad(xi,gpstruct,method,dx);
    hypw = 1;
else
    ymu = NaN(Nhyp,n);          % Observed value SD
    ys2 = NaN(Nhyp,n);           % Observed value SD
    fmu = NaN(Nhyp,n);          % Function mean
    fs2 = NaN(Nhyp,n);           % Function SD
    hypw = NaN(Nhyp,1);         % Sample weights
    
    if ~isfield(gpstruct,'x0') || isempty(gpstruct.x0)
        X = gpstruct.x;
    else
        xi = bsxfun(@minus,xi,gpstruct.x0);
        X = bsxfun(@minus,gpstruct.x,gpstruct.x0);
    end    
    
    % Compute prediction for each sample
    for iSample = 1:Nhyp        
        hyp = gpstruct.hyp(iSample);
        if isfield(gpstruct,'post') && ~isempty(gpstruct.post) && (numel(gpstruct.post) >= iSample)
            post = gpstruct.post(iSample);
        else
            post = [];
        end
        
        try
            [ymu(iSample,:),ys2(iSample,:),fmu(iSample,:),fs2(iSample,:)] = ...
                gppred1(X,xi,gpstruct,hyp,post);
                
            if isfield(gpstruct,'hypweight') && ~isempty(gpstruct.hypweight)
                hypw(iSample) = gpstruct.hypweight(iSample);
            else
                hypw(iSample) = 1;
            end
        catch
            % Failed sample is ignored
        end
    end
end

hypw = hypw./sum(isfinite(hypw(:)));    % Normalize samples

end

%GPPRED1 gp prediction for a single hyper-parameter sample
function [ymui,ys2i,fmui,fs2i] = gppred1(X,xi,gpstruct,hyp,post)

    if isfield(gpstruct,'marginalize') && gpstruct.marginalize ... % Marginalized GP
            && isfield(gpstruct,'hypHessian') && ~isempty(gpstruct.hypHessian) 
        % Pre-computed Laplace approximation at MAP
        %if isfield(gpstruct,'hypHessian') && ~isempty(gpstruct.hypHessian)
            hyp.Sigma_inv = gpstruct.hypHessian;
        %end
        
        if ~isempty(post)    
            [ymui,ys2i,fmui,fs2i] = mgp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,X,post,xi);
        else
            [ymui,ys2i,fmui,fs2i] = mgp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,X,gpstruct.y,xi);
        end
    else
        if ~isempty(post)
            [ymui,ys2i,fmui,fs2i] = mygp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,X,post,xi);
        else
            [ymui,ys2i,fmui,fs2i] = mygp(hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
        gpstruct.lik,X,gpstruct.y,xi);
        end
        
        % Recompute noiseless prediction if non-standard likelihood (e.g. warped gp)
        if numel(hyp.lik) > 1
            hyp.lik(end) = -Inf;
            [~, fmui, fs2i] = feval(gpstruct.lik{:}, hyp.lik, [], fmui, fs2i);
        end
    end
end