function uinit = initUniform(u0,LB,UB,PLB,PUB,Ninit,optimState,options)
%INITUNIFORM Initialize with points uniformly distributed in the box

nvars = numel(u0);

r = rand(nvars,Ninit);
uinit = bsxfun(@plus, PLB, bsxfun(@times, r, PUB-PLB));
