function [post nlZ dnlZ] = infGrid_Laplace(hyp, mean, cov, lik, x, y, opt)

% Laplace approximation to the posterior Gaussian process with covGrid
% covariance and (possibly) non-Gaussian likelihood.
% The (Kronecker) covariance matrix used is given by:
%   K = kron( kron(...,K{2}), K{1} ) = K_p x .. x K_2 x K_1.
%
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m and in conjunction with covGrid. See also infMethods.m.
%
% Copyright (c) by Hannes Nickisch 2015-03-05.
%
% See also INFMETHODS.M, COVGRID.M.

persistent last_alpha                                   % copy of the last alpha
if any(isnan(last_alpha)), last_alpha = zeros(size(last_alpha)); end   % prevent

inf = 'infLaplace';
cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
if ~strcmp(cov1,'covGrid'); error('Only covGrid supported.'), end    % check cov

xg = cov{3}; p = numel(xg);                                    % underlying grid
N = 1; D = 0; for i=1:p, N = N*size(xg{i},1); D = D+size(xg{i},2); end    % dims
[K,M,xe] = feval(cov{:}, hyp.cov, x);     % evaluate covariance mat constituents
for j=1:numel(K)                                               % expand Toeplitz
  if iscell(K{j}) && strcmp(K{j}{1},'toep'), K{j} = toeplitz(K{j}{2}); end
end
xe = M*xe; n = size(xe,1);
m = feval(mean{:}, hyp.mean, xe);                         % evaluate mean vector

likfun = @(f) feval(lik{:},hyp.lik,y,f,[],inf);        % log likelihood function

if nargin<=6, opt = []; end                        % make opt variable available
if isfield(opt,'cg_tol'),   cgtol = opt.cg_tol;       % stop conjugate gradients
else cgtol = 1e-5; end
if isfield(opt,'cg_maxit'), cgmit = opt.cg_maxit;      % number of cg iterations
else cgmit = min(n,20); end
if isfield(opt,'nlZ_exact'), exact = opt.nlZ_exact;        % use dense matrix ..
else exact = 0; end                             % .. computations instead of LCG

if any(size(last_alpha)~=[n,1])     % find a good starting point for alpha and f
  alpha = zeros(n,1);                      % start at mean if sizes do not match
else
  alpha = last_alpha;                                             % try last one
  if Psi(alpha,m,K,M,likfun) > -sum(likfun(m))   % default f==m better => use it
    alpha = zeros(n,1);
  end
end

% switch between optimisation methods
alpha = irls(alpha, m,K,M,likfun, opt);                       % run optimisation
alpha = irls(alpha, m,K,M,likfun, opt);                              % run again

f = M*kronmvm(K,M'*alpha) + m;                  % compute latent function values
last_alpha = alpha;                                     % remember for next call
[lp,dlp,d2lp,d3lp] = likfun(f); W = -d2lp;                 % evaluate likelihood
post.alpha = alpha;                            % return the posterior parameters
sW = sqrt(max(W,0)); post.sW = sW;             % remove sign in case of negative
mvm = @(t) mvmsWKsW(t,K,sW,M);                  % symmetric single parameter mvm
rep = @(sW,k) repmat(sW,1,size(k,2));
post.L = @(k) -rep(sW,k).*solveMVM(rep(sW,k).*k,mvm,cgtol,cgmit);

% diagnose optimality
err = @(x,y) norm(x-y)/max([norm(x),norm(y),1]);   % we need to have alpha = dlp
% dev = err(alpha,dlp);  if dev>1e-4, warning('Not at optimum %1.2e.',dev), end

if nargout>1
  if exact                               % exact marginal likelihood computation
    Ki = kronmvm(K,eye(N)); Ki = M*Ki*M';
    Li = chol(eye(n)+sW*sW'.*Ki);
    nlZ = alpha'*(f-m)/2 - sum(lp) + sum(log(diag(Li)));
  else                 % Fiedler's 1971 upper log determinant bound on A = I+K*W
    a = logdet_I_KW(K,W);
    nlZ = alpha'*(f-m)/2 - sum(lp) + a/2;                   % upper bound on nlZ
  end
end

if nargout>2                                           % do we want derivatives?
  dnlZ = hyp;                                   % allocate space for derivatives
  mvm = @(t) mvmsWKsW(t,K,sW,M);                          % single parameter mvm
  mvmZ = @(k) sW.*solveMVM(k.*sW,mvm,cgtol,cgmit);           % Z = inv(K+inv(W))
  if exact              % deriv. of ln|B| wrt W; g = diag(inv(inv(K)+diag(W)))/2
    A = Ki*diag(W)+eye(numel(m));
    Ci = Li'\(repmat(sW,1,numel(m)).*Ki);
    g = (diag(Ki)-sum(Ci.^2,1)')/2;
    dfhat = g.*d3lp;              % d nlZ / d fhat:  diag(inv(inv(K)+W)).*d3lp/2
  else
    g = zeros(n,1);
  end                                         % no implicit derivative if ~exact
  for i=1:length(hyp.cov)                                    % covariance hypers
    dK = feval(cov{:}, hyp.cov, x, [], i);
    for j=1:numel(dK)                                          % expand Toeplitz
      if iscell(dK{j})&&strcmp(dK{j}{1},'toep'), dK{j} = toeplitz(dK{j}{2}); end
    end
    if exact
      dKi = kronmvm(dK,eye(N)); dKi = M*dKi*M';
      dla = trace( A'\(diag(W)*dKi) );
    else
      h = 1e-4;                                % finite difference approximation
      hyp_cov_h = hyp.cov; hyp_cov_h(i) = hyp_cov_h(i)+h;
      K_h = feval(cov{:}, hyp_cov_h, x);
      for j=1:numel(K_h)                                       % expand Toeplitz
        if iscell(K_h{j})&&strcmp(K_h{j}{1},'toep')
          K_h{j} = toeplitz(K_h{j}{2});
        end
      end
      dla = (logdet_I_KW(K_h,W)-a)/h;          % finite difference approximation
    end
    dnlZ.cov(i) = dla/2 - alpha'*M*kronmvm(dK,M'*alpha)/2;       % explicit part
    if exact
      b = M*kronmvm(dK,M'*dlp);                        % inv(eye(n)+K*diag(W))*b
      dnlZ.cov(i) = dnlZ.cov(i)-dfhat'*(b-M*kronmvm(K,M'*mvmZ(b)));   % implicit
    end
  end
  for i=1:length(hyp.lik)                                    % likelihood hypers
    [lp_dhyp,dlp_dhyp,d2lp_dhyp] = feval(lik{:},hyp.lik,y,f,[],inf,i);
    dnlZ.lik(i) = -g'*d2lp_dhyp - sum(lp_dhyp);                  % explicit part
    if exact
      b = M*kronmvm(K,M'*(dlp_dhyp));                  % inv(eye(n)+K*diag(W))*b
      dnlZ.lik(i) = dnlZ.lik(i)-dfhat'*(b-M*kronmvm(K,M'*mvmZ(b)));   % implicit
    end
  end
  for i=1:length(hyp.mean)                                         % mean hypers
    dm = feval(mean{:}, hyp.mean, xe, i);
    dnlZ.mean(i) = -alpha'*dm;                                   % explicit part
    if exact
      dnlZ.mean(i) = dnlZ.mean(i)-dfhat'*(dm-M*kronmvm(K,M'*mvmZ(dm))); % implic
    end
  end
end

% Upper determinant bound on log |I+K*diag(W)| using Fiedler's 1971 inequality.
% K = kron( kron(...,K{2}), K{1} ), W = diag(w) both symmetric psd.
% The bound is exact for W = w*ones(N,1).
%
% Given nxn spd matrices C and D with ordered nx1 eigenvalues c, d 
% then exp(lb)=prod(c+d) <= det(C+D) <= prod(c+flipud(d))=exp(ub).
function [ub,lb] = logdet_I_KW(K,W)
  e = 1; for i=1:numel(K), e = kron(eigr(K{i}),e); end     % eigenvalue diagonal
  n = numel(W); N = numel(e);                                       % dimensions
  a = n/N*sort(e,'descend');                         % approximate spectrum of K
  if n>N, a = [a;zeros(n-N,1)]; else a = a(1:n); end
  ub = sum(log( 1 + a.*sort(W,'descend') ));     % See also Prob.III.6.14 in ..
  lb = sum(log( 1 + a.*sort(W, 'ascend') ));  % .. Matrix Analysis, Bhatia 1997.

% Evaluate criterion Psi(alpha) = alpha'*K*alpha + likfun(f), where 
% f = K*alpha+m, and likfun(f) = feval(lik{:},hyp.lik,y,  f,  [],inf).
function [psi,dpsi,f,alpha,dlp,W] = Psi(alpha,m,K,M,likfun)
  f = M*kronmvm(K,M'*alpha) + m;
  [lp,dlp,d2lp] = likfun(f); W = -d2lp;
  psi = alpha'*(f-m)/2 - sum(lp);
  if nargout>1, dpsi = M*kronmvm(K,M'*(alpha-dlp)); end

% Run IRLS Newton algorithm to optimise Psi(alpha).
function alpha = irls(alpha, m,K,M,likfun, opt)
  n = numel(alpha);
  if isfield(opt,'irls_maxit'), maxit = opt.irls_maxit; % max no of Newton steps
  else maxit = 20; end                                           % default value
  if isfield(opt,'irls_Wmin'),  Wmin = opt.irls_Wmin; % min likelihood curvature
  else Wmin = 0.0; end                                           % default value
  if isfield(opt,'irls_tol'),   tol = opt.irls_tol;     % stop Newton iterations
  else tol = 1e-6; end                                           % default value
  if isfield(opt,'cg_tol'),   cgtol = opt.cg_tol;     % stop conjugate gradients
  else cgtol = 1e-5; end
  if isfield(opt,'cg_maxit'), cgmit = opt.cg_maxit;    % number of cg iterations
  else cgmit = 2*n; end         
  
  smin_line = 0; smax_line = 2;           % min/max line search steps size range
  nmax_line = 10;                          % maximum number of line search steps
  thr_line = 1e-4;                                       % line search threshold
  Psi_line = @(s,alpha,dalpha) Psi(alpha+s*dalpha, m,K,M,likfun);  % line search
  pars_line = {smin_line,smax_line,nmax_line,thr_line};  % line seach parameters
  search_line = @(alpha,dalpha) brentmin(pars_line{:},Psi_line,5,alpha,dalpha);

  f = M*kronmvm(K,M'*alpha) + m;
  [lp,dlp,d2lp] = likfun(f); W = -d2lp;
  Psi_new = Psi(alpha,m,K,M,likfun);
  Psi_old = Inf;  % make sure while loop starts by the largest old objective val
  it = 0;                          % this happens for the Student's t likelihood
  while Psi_old - Psi_new > tol && it<maxit                       % begin Newton
    Psi_old = Psi_new; it = it+1;                               % limit stepsize
    W = max(W,Wmin); sW = sqrt(W);     % reduce steps by incr curvature of bad W
    b = W.*(f-m) + dlp; mvm = @(t) mvmsWKsW(t,K,sW,M);    % single parameter mvm
    dalpha = sW.*solveMVM(b./sW,mvm,cgtol,cgmit) - alpha; % Newt dir+line search
%     c = sW.*(M*kronmvm(K,M'*b)); % like in GPML book
%     dalpha = b - sW.*solveMVM(c,mvm,cgtol,cgmit) - alpha;
    [s_line,Psi_new,n_line,dPsi_new,f,alpha,dlp,W] = search_line(alpha,dalpha);
  end                                                  % end Newton's iterations

% solve q = mvm(p) via conjugate gradients
function q = solveMVM(p,mvm,tol,maxit)
  [q,flag,relres,iter] = conjgrad(mvm,p,tol,maxit);                   % like pcg
% if ~flag,error('Not converged after %d iterations, r=%1.2e\n',iter,relres),end
%   maxres = 0;
%   for i=1:size(p,2)                             % TODO: maybe have a better CG
%     [txt,q(idx,i),flag,relres] = evalc('pcg(mvm,p(idx,i),tol,maxit)');
%     maxres = max(maxres,relres);
%   end
%   if maxres>tol, error('Not converged, r=%1.2e\n',maxres),end

% symmetric mvm so that q = diag(sW)*M*K*M'*diag(sW)*p + p
% using Kronecker representation
function q = mvmsWKsW(p,K,sW,M)
  q = repmat(sW,1,size(p,2)).*p;
  q = M*kronmvm(K,M'*q);
  q = repmat(sW,1,size(p,2)).*q + p;

%------------------------------------------------------------------------------%

% Solve x=A*b with symmetric A(n,n), b(n,m), x(n,m) using conjugate gradients.
% The method is along the lines of PCG but suited for matrix inputs b.
function [x,flag,relres,iter,r] = conjgrad(A,b,tol,maxit)
if nargin<3, tol = 1e-10; end
if nargin<4, maxit = min(size(b,1),20); end
x0 = zeros(size(b)); x = x0;
if isnumeric(A), r = b-A*x; else r = b-A(x); end, r2 = sum(r.*r,1); r2new = r2;
nb = sqrt(sum(b.*b,1)); flag = 0; iter = 1;
relres = sqrt(r2)./nb; todo = relres>=tol; if ~any(todo), flag = 1; return, end
on = ones(size(b,1),1); r = r(:,todo); d = r;
for iter = 2:maxit
  if isnumeric(A), z = A*d; else z = A(d); end
  a = r2(todo)./sum(d.*z,1);
  a = on*a;
  x(:,todo) = x(:,todo) + a.*d;
  r = r - a.*z;
  r2new(todo) = sum(r.*r,1);
  relres = sqrt(r2new)./nb; cnv = relres(todo)<tol; todo = relres>=tol;
  d = d(:,~cnv); r = r(:,~cnv);                           % get rid of converged
  if ~any(todo), flag = 1; return, end
  b = r2new./r2;                                               % Fletcher-Reeves
  d = r + (on*b(todo)).*d;
  r2 = r2new;
end

% Perform a matrix vector multiplication b = A*x with a matrix A being a
% Kronecker product given by A = kron( kron(...,As{2}), As{1} ).
function b = kronmvm(As,x,transp)
if nargin>2 && ~isempty(transp) && transp   % transposition by transposing parts
  for i=1:numel(As), As{i} = As{i}'; end
end
m = zeros(numel(As),1); n = zeros(numel(As),1);                  % extract sizes
for i=1:numel(n)
  if iscell(As{i}) && strcmp(As{i}{1},'toep')
    m(i) = size(As{i}{2},1); n(i) = size(As{i}{2},1);
  else [m(i),n(i)] = size(As{i});
  end
end
d = size(x,2);
b = x;
for i=1:numel(n)
  a = reshape(b,[prod(m(1:i-1)), n(i), prod(n(i+1:end))*d]);    % prepare  input
  if iscell(As{i}) && strcmp(As{i}{1},'toep')        % apply along one dimension
    b = zeros([prod(m(1:i-1)), m(i), prod(n(i+1:end))*d]);      % prepare output
    for j=1:size(b,3), b(:,:,j) = toepmvmsym(As{i}{2},a(:,:,j)',true)'; end
  else
    tmp = reshape(permute(a,[1,3,2]),[],n(i))*As{i}';
    b = permute(reshape(tmp,[size(a,1),size(a,3),m(i)]),[1,3,2]);
  end
end
b = reshape(b,prod(m),d);                        % bring result in correct shape

% Real eigenvalues and eigenvectors up to the rank of a real symmetric matrix.
% Decompose A into V*D*V' with orthonormal matrix V and diagonal matrix D.
% Diagonal entries of D obave the rank r of the matrix A as returned by
% the call rank(A,tol) are zero.
function [V,D] = eigr(A,tol)
[V,D] = eig((A+A')/2); n = size(A,1);    % decomposition of strictly symmetric A
d = max(real(diag(D)),0); [d,ord] = sort(d,'descend');        % tidy up and sort
if nargin<2, tol = size(A,1)*eps(max(d)); end, r = sum(d>tol);     % get rank(A)
d(r+1:n) = 0;                             % set junk eigenvalues to strict zeros
if nargout==1, V = d; return, end
D = diag(d);                                      % inflate to respect interface
V(:,1:r) = real(V(:,ord(1:r))); V(:,r+1:n) = null(V(:,1:r)'); % ortho completion