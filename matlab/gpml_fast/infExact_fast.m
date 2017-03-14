function [post nlZ dnlZ] = infExact_fast(hyp, mean, cov, lik, x, y)

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
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood');
end
 
[n, D] = size(x);
if nargout > 2                                      % do we want derivatives?
  [K,dK] = feval(cov{:}, hyp.cov, x);                 % evaluate covariance matrix and derivatives
else
  K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
end
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss

smallnoise = sn2 < 1e-6;       % very tiny sn2 can lead to numerical trouble

if smallnoise
    M = K+sn2*eye(n);           % Covariance with noise
else
    M = K/sn2+eye(n);           % B matrix
end
[L,p] = chol(M);                % Try computing Cholesky factor

% M
%[V,D] = eig(M);
% V(:,1)'
%if D(1,1) < 0
%    M/max(M(:))
%   diag(D) 
%   V(:,1)
%   sn2
%   pause
%end

%[D(1,1)/D(2,2)]
%[D(1,1),D(2,2),sn2]

% log10(diag(lambda))'
% max(lambda)/min(lambda)

if p~=0     % Failed Cholesky decomposition, compute nearest SPD matrix
    error('a');
    Mold = M;
    
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
        if k > 5; error('Cannot compute Cholesky decomposition.'); end        
        lambda = eig(M);
        
        mineig = min(lambda);

        mineigold = log10(abs(min(eig(Mold))));
       
        [log10(abs(max(lambda))),log10(abs(mineig)),mineigold,log10(sn2)]
        kappa = 0.125;       % Sometimes even a small nudge is sufficient
        % kappa = 0.125;      % Sometimes even a small nudge is sufficient
        % M = M + (-kappa*mineig*k.^2 + eps(mineig))*eye(size(M));        
        % M = M + (-kappa*mineig*k.^2 + sqrt(eps)*mineig)*eye(size(M));        
        M = M + kappa*abs(mineig)*k.^2*eye(size(M));
        % M = M + kappa*mineig*k.^2*eye(size(M));
        % k, [lambda';eig(M)']
        [L,p] = chol(M);
    end
    
    % [sort(eig(M))'; sort(eig(Mold))']
    
end

if smallnoise
    sl = 1;
    pL = -solve_chol(L,eye(n));     % L = -inv(K+inv(sW^2))
else
    sl = sn2;
    pL = L;                         % L = chol(eye(n)+sW*sW'.*K)
end

alpha = solve_chol(L,y-m)/sl;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                  % sqrt of noise precision vector
post.L = pL;

if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sl)/2;   % -log marg lik
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;                                 % allocate space for derivatives
    Q = solve_chol(L,eye(n))/sl - alpha*alpha';     % precompute for convenience
    for i = 1:numel(hyp.cov)
      dnlZ.cov(i) = sum(sum(Q.*dK(:,:,i)))/2;
    end
    dnlZ.lik = sn2*trace(Q);
    for i = 1:numel(hyp.mean)
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*alpha;
    end
  end
end
