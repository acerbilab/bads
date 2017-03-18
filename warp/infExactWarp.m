function [post nlZ dnlZ] = infExactWarp(hyp, mean, cov, lik, x, y)

% Exact inference for a GP with Gaussian likelihood. Compute a parametrization
% of the posterior, the negative log marginal likelihood and its derivatives
% w.r.t. the hyperparameters. See also "help infMethods".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2015-07-13.
%                                      File automatically generated using noweb.
%
% See also INFMETHODS.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGaussWarpExact')               % NOTE: no explicit call to likGaussWarpExact
  error('Exact inference only possible with warped exact Gaussian likelihood');
end
warp = lik{2}; if ischar(warp); warp = str2func(warp); end

ng = feval(warp{:});               % number of hyperparameters for the warping function
nhyp = feval(lik{:});              % number of hyperparameters
% if nargin<4, varargout = {nhyp}; return, end       % report number of parameters
nhyp = eval(nhyp);
if nhyp>length(hyp.lik), error('not enough hyperparameters'), end

[gy,lgpy] = feval(warp{:},y,hyp.lik(1:ng));                     % evaluate warping function

% figure(100); scatter(y,gy); drawnow;

[n, D] = size(x);
K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
m = feval(mean{:}, hyp.mean, x);                    % evaluate mean vector
[gm,lgpm] = feval(warp{:},m,hyp.lik(1:ng));         % evaluate warped mean

sn2 = exp(2*hyp.lik(end));                          % noise variance of likGauss
if sn2<1e-6                        % very tiny sn2 can lead to numerical trouble
  L = chol(K+sn2*eye(n)); sl =   1;   % Cholesky factor of covariance with noise
  pL = -solve_chol(L,eye(n));                            % L = -inv(K+inv(sW^2))
else
  L = chol(K/sn2+eye(n)); sl = sn2;                       % Cholesky factor of B
  pL = L;                                           % L = chol(eye(n)+sW*sW'.*K)
end
alpha = solve_chol(L,gy-gm)/sl;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                  % sqrt of noise precision vector
post.L = pL;

if nargout>1                               % do we want the marginal likelihood?
  nlZ = (gy-gm)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sl)/2 - sum(lgpy);   % -log marg lik
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;                                 % allocate space for derivatives
    Q = solve_chol(L,eye(n))/sl - alpha*alpha';     % precompute for convenience
    for i = 1:numel(hyp.cov)
      dnlZ.cov(i) = sum(sum(Q.*feval(cov{:}, hyp.cov, x, [], i)))/2;
    end
    for i = 1:ng
        [dgy,dlgpy] = feval(warp{:},y,hyp.lik(1:ng),i);
        dgm = feval(warp{:},m,hyp.lik(1:ng),i);
        dnlZ.lik(i) = -sum(dlgpy) + (dgy - dgm)'*alpha;
    end
    dnlZ.lik(ng+1) = sn2*trace(Q);
    for i = 1:numel(hyp.mean)
      dnlZ.mean(i) = -(exp(lgpm)'.*feval(mean{:}, hyp.mean, x, i)')*alpha;
    end
  end
  
end