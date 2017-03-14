function xs = searchGrid(x,gpstruct,LB,UB,optimState,options)
%SEARCHGRID Full grid search.

if nargin < 1
    xs = 'grd';
    return;
end

nvars = numel(x);
iter = optimState.iter;
MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

try
    r = i4_sobol_generate(nvars,options.Nsearch,(iter-1)*options.Nsearch)';
catch
    r = rand(options.Nsearch,nvars);
end

r = 2*(r - 0.5);

%landscapeplot(@(x_) fun(x_.*optimState.scale + optimState.x0),x,LB,UB,SearchFactor*MeshSize,gpstruct,v,31);
%drawnow;

xs = bsxfun(@plus,x,bsxfun(@times,r,MeshSize*SearchFactor));

end