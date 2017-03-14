function u = force2grid(u,optimState,tol)
%FORCE2GRID 
if nargin < 3 || isempty(tol); tol = optimState.searchmeshsize; end

u = tol.*round(u./tol);

end