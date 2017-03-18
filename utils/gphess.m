function [Hfmu] = gphess(x0,gpstruct,dx,method)
%GPHESS Hessian of Gaussian process mean prediction

if nargin < 3 || isempty(dx); dx = 1e-6; end
if nargin < 4 || isempty(method); method = 'central'; end

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
    
    switch lower(method)
        case 'five-points'
            xi = bsxfun(@plus,x0,dx/12*[2*eye(D);eye(D);-eye(D);-2*eye(D)]);
    end
else
    xi = [];
end

xi = [xi; x0];

% Compute predictions and derivatives for each hyperparameter sample
for i = 1:Nhyp
    
    try
        [ymui,ys2i,fmui,fs2i] = gp(gpstruct.hyp(i),gpstruct.inf,gpstruct.mean,gpstruct.cov,gpstruct.lik,gpstruct.x,gpstruct.y,xi);
    catch
        ymui = NaN(size(xi,1),1); ys2i = NaN(size(xi,1),1);
        fmui = NaN(size(xi,1),1); fs2i = NaN(size(xi,1),1);
    end

    % Compute prediction
    ymu(i) = ymui(end);
    fmu(i) = fmui(end);
    ys2(i) = ys2i(end);
    fs2(i) = fs2i(end);
    
    % Compute derivatives
    if nargout > 4
        switch method
            case 'central'
                dymu(i,:) = (ymui(1:D) - ymui(D+(1:D)))./dx;
                dfmu(i,:) = (fmui(1:D) - fmui(D+(1:D)))./dx;
                dys2(i,:) = (ys2i(1:D) - ys2i(D+(1:D)))./dx;
                dfs2(i,:) = (fs2i(1:D) - fs2i(D+(1:D)))./dx;
            case 'forward'
                dymu(i,:) = (ymui(1:D) - ymui(end))./dx;
                dfmu(i,:) = (fmui(1:D) - fmui(end))./dx;
                dys2(i,:) = (ys2i(1:D) - ys2i(end))./dx;
                dfs2(i,:) = (fs2i(1:D) - fs2i(end))./dx;                
            case 'backward'
                dymu(i,:) = -(ymui(1:D) - ymui(end))./dx;
                dfmu(i,:) = -(fmui(1:D) - fmui(end))./dx;
                dys2(i,:) = -(ys2i(1:D) - ys2i(end))./dx;
                dfs2(i,:) = -(fs2i(1:D) - fs2i(end))./dx;                
            case 'five-points'
                dymu(i,:) = (-ymui(1:D) + 8*ymui(D+(1:D)) - 8*ymui(2*D+(1:D)) + ymui(3*D+(1:D)))./dx;
                dfmu(i,:) = (-fmui(1:D) + 8*fmui(D+(1:D)) - 8*fmui(2*D+(1:D)) + fmui(3*D+(1:D)))./dx;
                dys2(i,:) = (-ys2i(1:D) + 8*ys2i(D+(1:D)) - 8*ys2i(2*D+(1:D)) + ys2i(3*D+(1:D)))./dx;
                dfs2(i,:) = (-fs2i(1:D) + 8*fs2i(D+(1:D)) - 8*fs2i(2*D+(1:D)) + fs2i(3*D+(1:D)))./dx;
                
        end
    end   
end

end