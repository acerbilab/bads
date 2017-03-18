function Bnew = pollGPS2N(B,x,gpstruct,LB,UB,optimState,options)
%POLLGPS2N Poll 2N fixed basis vectors (Generalized Pattern Search).

nvars = numel(x);

if isempty(B)
    % Axis-oriented poll vector set
    Bnew = [eye(nvars); -eye(nvars)];
else
    Bnew = [];
end

end