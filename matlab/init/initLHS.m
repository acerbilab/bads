function [uinit,bestscore] = initLHS(u0,LB,UB,PLB,PUB,Ninit,optimState,options)
%INITLHS Initialize with Latin hypercube sample design in the box

nvars = size(u0,2);

% Only consider periodic variables that cover full range
per = optimState.periodicvars;
for d = find(per)
    if PLB(d) > LB(d) || PUB(d) < UB(d); per(d) = 0; end
end

maxiter = 1e3;

% Compute LHS sample
[uinit,bestscore] = lhs(Ninit,nvars,PLB,PUB,u0,maxiter,per);

end