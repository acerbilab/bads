function uinit = initSobol(u0,LB,UB,PLB,PUB,Ninit,optimState,options)
%INITSOBOL Initialize with space-filling design in the box

nvars = numel(u0);

try
    % Try to generate a quasirandom Sobol space-filling sequence
    
    MaxSeed = 997;    
    if all(isfinite(u0))                        % Seed depends on U0
        strseed = num2str(u0(1:min(10,end)));       
        seed = mod(prod(uint64(strseed)),MaxSeed)+1;
    else
        seed = randi(MaxSeed);                  % Random seed
    end
    r = i4_sobol_generate(nvars,Ninit,seed)';
    uinit = bsxfun(@plus, PLB, bsxfun(@times, r, PUB-PLB));
catch
    % Latin hypercube sampling if generation of Sobol sequence failed
    uinit = initLHS(u0,LB,UB,PLB,PUB,Ninit,optimState,options);
end