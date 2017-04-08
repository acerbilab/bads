function [X,bestscore] = lhs(N,p,LB,UB,x0,maxiter,per)
%LHS Latin hypercube sample
%
%  X = LHS(N,P) creates a Latin hypercube sample design with N points on 
%  the unit hypercube [0,1]^P.
%
%  X = LHS(N,P,LB,UB) creates a hypercube with lower/upper bounds LB and UB.
%
%  X = LHS(N,P,LB,UB,X0) assumes that the cell with X0 is already occupied.
%
%  X = LHS(N,P,LB,UB,X0,MAXITER) optimizes the lhs design up to MAXITER 
%  times (default MAXITER=100).
%
% Version 0.8
% Author:       Luigi Acerbi
% Release Date: Jan/22/2016

if nargin < 3; LB = []; end
if nargin < 4; UB = []; end
if nargin < 5; x0 = []; end
if nargin < 6 || isempty(maxiter); maxiter = 100; end
if nargin < 7 || isempty(per); per = zeros(1,p); end

if isempty(LB); LB = zeros(1,p); end
if isempty(UB); UB = ones(1,p); end

bestscore = -Inf;

if ~isempty(x0)
    % Assign initial point to LHS grid
    x0 = (x0 - LB)./(UB-LB);
    i0 = floor((N+1)*x0);
    i0 = max(min(i0, N), 0);            % Enforce bounds
else
    i0 = [];
end

% Compute period
f = find(per);
per(f) = UB(f) - LB(f);         % Finite (periodic)
per(~f) = Inf;                  % Infinite (non-periodic)

% Compute several grids until iterations run out
while maxiter > 0
    % Compute locally optimized grid
    [Xnew,b,iter] = producegrid(N,p,maxiter,x0,i0,per);
    
    % Compare to previous grid, substitute if better
    if b > bestscore
        X = Xnew;
        bestscore = b;
    end
    
    % Decrease available iterations
    maxiter = maxiter - iter;
end

% Rescale grid
X = bsxfun(@plus, LB, bsxfun(@times, X, UB-LB));

end

%--------------------------------------------------------------------------
function [X,bestscore,iter] = producegrid(N,p,maxiter,x0,i0,per)
%PRODUCEGRID Produce a local LHS grid

% Start with a plain LHS sample over a grid
X = getsample(N,p,i0);
iter = 1;
bestscore = mindist([x0;X],per);

if N == 1; return; end  % Only one point, return

lastimprovement = 0;
n1 = 1 + randi(N-1);    % First input

% Try local improvement by swapping array elements
while iter < maxiter
    iter = iter + 1;
    
    col = mod(iter,p)+1;    % Choose variable
    
    % n1 = 1 + randi(Ninit-1);    % First input (do not change)
    if N > 2
        n2 = 1+randi(N-2);      % Randomly chosen second input
    else
        n2 = 1;
    end
    if n2 >= n1; n2 = n2 + 1; end
    Xnew = X;
    temp = Xnew(n1,col);
    Xnew(n1,col) = Xnew(n2,col);
    Xnew(n2,col) = temp;

    [newscore,n1,n2] = mindist([x0;Xnew],per);
    if ~isempty(x0)
        if n1 == 1; n1 = n2-1; else n1 = n1-1; end
    end        
    if newscore > bestscore
        X = Xnew;
        bestscore = newscore;
        lastimprovement = iter;
        n1 = 1 + randi(N-1);
        % [iter bestscore]
    end

    if iter - lastimprovement > 3*p
        n1 = 1 + randi(N-1);    % Randomize first input
    end
    
    if iter - lastimprovement > 50
        % Quit if more than fifty iterations with no improvement
        break;
    end
end

end

%--------------------------------------------------------------------------
function x = getsample(n,p,i0)
%GETSAMPLE Get a new LHS grid

if nargin < 3; i0 = []; end
% if isempty(i0); n = n-1; end

x = rand(n,p);
for i=1:p
    [~, rowidx] = sort(x(:,i));
    r(rowidx) = 0:n-1;
    x(:,i) = r(:);                  % rank x
end

nf = n + 1;
if ~isempty(i0)
    x = x + bsxfun(@ge, x, i0);
end

x = x + rand(size(x));
x = x / nf;

end


%--------------------------------------------------------------------------
function [m2,row,col] = mindist(X,per)
%MINDIST Compute minimum distance between vectors in X

% Compute distance (takes into account periodic variables)
x1(:,:,1) = X;
x2(1,:,:) = X';
d2(:,:) = dist(x1,x2,per);

d2 = tril(d2,-1);
d2(d2 == 0) = Inf;
[a,row] = min(d2);
[m2,col] = min(a);
row = row(col);

% m2 = min(d2(d2 > 0));

end

%--------------------------------------------------------------------------
function d2 = dist(x1,x2,per)
%DIST Vector distance

f = isfinite(per);

if any(f)
    D = bsxfun(@minus, x1, x2);
    for d = find(f)
        D(:,d,:) = min(abs(D(:,d,:)), per(d) - abs(D(:,d,:)));
    end    
else
    D = bsxfun(@minus, x1, x2);
end

d2 = sum(D.^2, 2);

end