function x = periodCheck(x,LB,UB,optimState,xgridscaled)
%PERIODCHECK Check and move periodic variables to main range.

if ~any(optimState.periodicvars); return; end

if nargin < 5 || isempty(xgridscaled); xgridscaled = 0; end

% If X is already rescaled, rescale also LB and UB
if xgridscaled
    UB = gridunits(UB,optimState);
    LB = gridunits(LB,optimState);
end

period = UB - LB;

for d = find(optimState.periodicvars)
    x(:,d) = mod(x(:, d) - LB(d), period(d)) + LB(d);    
end

end
