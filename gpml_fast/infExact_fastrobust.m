function [post,nlZ,dnlZ] = infExact_fastrobust(kmax, hyp, mean, cov, lik, x, y, s)

% Fast, robust inference for a GP with Gaussian likelihood. Compute a 
% parametrization of the posterior, the negative log marginal likelihood 
% and its derivatives w.r.t. the hyperparameters. See also "help infMethods".
% 
% To be used in combination with 'fast' covariance functions (e.g.
% covSEard_fast and such). 
%
% This implementation also handles badly conditioned covariance matrices, 
% which sometimes become non-positive definite and cause a failure of the 
% Cholesky decomposition. In these cases, the covariance matrix is corrected 
% to the nearest symmetric positive semidefinite (SPD) matrix in Frobenius 
% norm (see [1]; based on an implementation by John d'Errico).
%
% [1] Higham NJ, "Computing a nearest symmetric positive semidefinite 
% matrix", Linear Algebra Appl, 1988.
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% Added support for 'fast' inference and robust Cholesky decomposition by
% Luigi Acerbi, 2016-01-03.
% Original code by Carl Edward Rasmussen and Hannes Nickisch, 2015-07-13.
%
% See also INFMETHODS.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if strcmp(likstr,'likGauss')    % NOTE: no explicit call to likGauss
    henoise_flag = false;
elseif strcmp(likstr,'likGaussHe')
    henoise_flag = true;        % Heteroskedastic noise    
else
    error('Exact inference only possible with Gaussian likelihood');
end

if isempty(kmax); kmax = 5; end
 
[n, D] = size(x);
if nargout > 2                                      % do we want derivatives?
  [K,dK] = feval(cov{:}, hyp.cov, x);                 % evaluate covariance matrix and derivatives
else
  K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
end
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector

sn2_base = exp(2*hyp.lik);                               % noise variance of likGauss

if henoise_flag && ~isempty(s)
    % if isempty(s); error('Input-dependent vector S is empty.'); end
    sn2 = sn2_base + s.^2;                               % Vector of observation variance
else
    sn2 = sn2_base;
end

Lchol = min(sn2) >= 1e-6;       % tiny sn2 can lead to numerical trouble

% smallnoise = sn2 < 1e-6;       
if Lchol
    if isscalar(sn2)
        sn2div = sn2;
        sn2_mat = eye(n);
    else
        sn2div = min(sn2);
        sn2_mat = diag(sn2/sn2div);
    end    
    
    M = K/sn2div+sn2_mat;       % B matrix
else
    if isscalar(sn2)
        sn2_mat = sn2*eye(n);
    else
        sn2_mat = diag(sn2);
    end
    M = K+sn2_mat;           % Covariance with noise
end

[L,p] = chol(M);                % Try computing Cholesky factor

if p~=0     % Failed Cholesky decomposition, compute nearest SPD matrix
    if kmax <= 0; error('Cannot compute Cholesky decomposition.'); end
    % Mold = M;
    
    M = (M + M')/2;         % Ensure M is symmetric
    [U,Sigma,V] = svd(M);
    H = V*Sigma*V';         % Symmetric polar factor H is SPD
    M = (M+H)/2;            
    M = (M + M')/2;         % Ensure symmetry again
    [L,p] = chol(M);        % Retry Cholesky decomposition
    k = 0;
    while p ~= 0            % Failed again, add a small diagonal component
        k = k + 1;
        % We do not want to change the matrix too much
        if k > kmax-1; error('Cannot compute Cholesky decomposition.'); end        
        lambda = eig(M);        
        mineig = min(lambda);       
        kappa = 0.05;        % Sometimes even a small nudge is sufficient
        M = M + kappa*abs(mineig)*k.^2*eye(size(M));
        [L,p] = chol(M);
    end
    
    % [sort(eig(M))'; sort(eig(Mold))']
    
end

if Lchol
    sl = sn2div;
    pL = L;                         % L = chol(eye(n)+sW*sW'.*K)
else
    sl = 1;
    pL = -solve_chol(L,eye(n));     % L = -inv(K+inv(sW^2))
end

alpha = solve_chol(L,y-m)/sl;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(min(sn2));            % sqrt of noise precision vector
post.L = pL;
post.Lchol = Lchol;

if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sl)/2;   % -log marg lik
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;                                 % allocate space for derivatives
    Q = solve_chol(L,eye(n))/sl - alpha*alpha';     % precompute for convenience
    for i = 1:numel(hyp.cov)
      dnlZ.cov(i) = sum(sum(Q.*dK(:,:,i)))/2;
    end
    dnlZ.lik = sn2_base*trace(Q);
    for i = 1:numel(hyp.mean)
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*alpha;
    end
  end
end
