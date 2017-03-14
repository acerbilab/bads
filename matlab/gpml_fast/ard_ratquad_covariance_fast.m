% ARD_RATQUAD_COVARIANCE rational quadratic covariance with ARD.
%
% This provides a GPML-compatible covariance function implementing the
% rational quadratic covariance with automatic relevance
% determination (ARD). This can be used as a drop-in replacement for
% covRQard.
%
% This implementation supports an extended GPML syntax that allows
% calculating the Hessian of K with respect to any pair of
% hyperparameters. The syntax is:
%
%   dK2_didj = ard_ratquad_covariance(theta, x, z, i, j)
%
% where dK2_didj is \partial^2 K / \partial \theta_i \partial \theta_j,
% and the Hessian is evalauted at K(x, z). As in the derivative API,
% if z is empty, then the Hessian is evaluated at K(x, x).  Note that
% the option of setting z = 'diag' for Hessian computations is not
% supported due to no obvious need.
%
% These Hessians can be used to ultimately compute the Hessian of the
% GP training likelihood (see, for example, exact_inference.m).
%
% The hyperparameters are the same as for covRQard.
%
% See also COVRQARD, COVFUNCTIONS.

% Copyright (c) 2015 Luigi Acerbi.

function [K,dK,HK] = ard_ratquad_covariance_fast(theta, x, z, i, j)

  % used during gradient and Hessian calculations to avoid constant recomputation
  % persistent D2;
  % persistent K;

  % call covRQard for everything but Hessian calculation
  if (nargin <= 1)
    K = covRQard_fast;
  elseif (nargout == 2) && nargin < 4
      if (nargin == 2)
        [K,dK] = covRQard_fast(theta, x);
      elseif (nargin == 3)
        [K,dK] = covRQard_fast(theta, x, z);
      end
  elseif nargout == 1
      if (nargin == 2)
        K = covRQard_fast(theta, x);
      elseif (nargin == 3)
        K = covRQard_fast(theta, x, z);
      elseif (nargin == 4)
        K = covRQard_fast(theta, x, z, i);
      elseif (nargin == 5)
        K = ard_ratquad_covariance(theta, x, z, i, j);
      end
      
  % Hessian with respect to \theta_i \theta_j
  elseif nargout == 3

    % ensure i <= j by exploiting symmetry
    %if (i > j)
    %  result = ard_ratquad_covariance(theta, x, z, j, i);
    %  return;
    %end
    
    % needed to compute number of hyperparameters
    D = size(x, 2);

    [K,dK] = covRQard_fast(theta, x, z);

    xeqz = isempty(z);
    ell = exp(-theta(1:D));
    sf2 = exp(2*theta(D+1));
    alpha = exp(theta(D+2));
    
    if xeqz                                       % symmetric matrix Kxx
      D2 = sq_dist(diag(ell)*x');
    else                                          % cross covariances Kxz
      D2 = sq_dist(diag(ell)*x',diag(ell)*z');
    end        
    M = (1+0.5*D2/alpha);
    
    C = sq_dist_fast(x,[],ell);
    HK = zeros(size(x,1),size(x,1),D+2,D+2);

    % Hessian involving length scales only
    T = sf2*M.^(-alpha-2).*((alpha+1)/alpha);
    for j = 1:D
        for i = 1:j-1
            HK(:,:,i,j) = T .* C(:,:,i) .* C(:,:,j);
        end
    end
    
    T = sf2*bsxfun(@minus, ...
        bsxfun(@times, M.^(-alpha-2).*((alpha+1)/alpha), C), ...
        2*M.^(-alpha-1));
    for i = 1:D
        HK(:,:,i,i) = T(:,:,i).*C(:,:,i);
    end
                
    % Hessians involving the log output scale
    j = D+1;
    for i = 1:D+1
        HK(:,:,i,j) = 2*dK(:,:,i);
    end
    i = D+1;
    j = D+2;
    HK(:,:,i,j) = 2*dK(:,:,j);
        
    % Hessian involving length scale and shape parameter alpha
    j = D+2;
    HK(:,:,1:D,j) = bsxfun(@times, ...
        C, ...
        sf2*M.^(-alpha-1).*(-log(M)*alpha + (alpha+1)*(1 - 1./M)));
        
    % Hessian involving shape parameter alpha only
    i = D+2;
    j = D+2;        
    HK(:,:,i,j) = sf2*alpha./((D2 + 2*alpha).^2) .* M.^(-alpha) .* ...
        (D2.*(2*alpha + D2*(2+alpha)) + ...
        (D2 + 2*alpha).*log(M) .* ...
        (-D2 -2*(1+D2)*alpha + alpha*(D2 + 2*alpha).*log(M)));
    
    HK(isnan(HK)) = 0;    
  end