function [xs,optimState] = searchNewton(x,gpstruct,LB,UB,optimState,options)
%SEARCHNEWTON Newton search along GP derivative.

if nargin < 1
    xs = 'lin';
    return;
end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

% Compute GP derivative at x
% [~,~,~,~,~,~,dfmu,~] = gpgrad(x,gpstruct,'central',Step);
% if any(isnan(dfmu)); dfmu = randn(size(x)); end

% Compute Newton vector
v = (-optimState.Binv*optimState.grad)';    
v = v./sqrt(v*v');                      % Normalize vector

if 0
    r = 2*(rand(options.Nsearch,1) - 0.5);
else
    r = 2 + 2*randn(options.Nsearch,1);
    % warning('Cambiato searchNewton.')
end
alpha = SearchFactor*MeshSize*r;

%fun = optimState.fun;
%landscapeplot(@(x_) fun(x_.*optimState.scale + optimState.x0),x,LB,UB,SearchFactor*MeshSize,gpstruct,v,31);
%drawnow;

% Points along the line
xs = bsxfun(@plus,x,bsxfun(@times,alpha,v));

end