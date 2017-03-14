% EXACT_INFERENCE infExact replacement supporting Hessian calculation.
%
% This provides a GPML-compatible inference method performing exact
% inference with a Gaussian likelihood that supports an extended API
% allowing the calculation of:
%
% - the partial derivatives of \alpha with respect to \theta,
% - the partial derivatives of diag W^{-1} with respect to \theta, and
% - the Hessian of the negative log likelihood at \theta.
%
% The former two can be used to compute the gradient of the latent
% predictive mean and variance of the approximate posterior GP with
% respect to the hyperparameters.
%
% This can be used as a drop-in replacement for infExact with no extra
% computational cost.
%
% Usage
% -----
%
% The API is identical to GPML inference methods, expect for three
% additional optional arguments:
%
%   [posterior, nlZ, dnlZ, dalpha, dWinv, HnlZ] = ...
%       exact_inference(theta, mean_function, covariance_function, ...
%                       likelihood, x, y);
%
% dalpha and dWinv provide the partial derivatives of the posterior
% parameters \alpha and W^{-1} with respect to \theta. These
% arrangement of these structs is similar to the dnlZ struct. For
% example, dalpha.cov(:, 1) gives the derivative of \alpha with
% respect to the first covariance hyperparameter.
%
% HnlZ provides the Hessian of -\log p(y | X, \theta). See hessians.m
% for information regarding the Hessian struct HnlZ.
%
% Requirements
% ------------
%
% The posterior derivatives dalpha and dWinv can be used with
% unmodified GPML mean, covariance, and likelihood functions.
%
% To compute the Hessian HnlZ, both the mean and covariance functions
% must support an extended GPML syntax that allows for calculating the
% Hessian of the training mean mu or training covariance K with
% respect to any pair of hyperparameters. The syntax for mean
% functions is:
%
%   d2mu_didj = mean_function(theta, x, i, j);
%
% where d2mu_didj is \partial^2 mu(x) / \partial \theta_i \partial \theta_j.
%
% The syntax for covariance functions is similar:
%
%   d2K_didj = covariance_function(theta, x, [], i, j);
%
% where d2K_didj is \partial^2 K(x, x) / \partial \theta_i \partial \theta_j.
%
% See also INFMETHODS, HESSIANS.

% Copyright (c) 2013--2015 Roman Garnett.

function [posterior, nlZ, dnlZ, dalpha, dWinv, HnlZ] = ...
      exact_inference_robust(theta, mean_function, covariance_function, ...
                      ~, x, y)

  % If additional outputs are not requested, simply call infExact and
  % return.
  if (nargout <= 1)
    posterior = ...
        infExact_fast(theta, mean_function, covariance_function, 'likGauss', x, y);
    return;
  elseif (nargout == 2)
    [posterior, nlZ] = ...
        infExact_fast(theta, mean_function, covariance_function, 'likGauss', x, y);
    return;
  elseif (nargout == 3)
    [posterior, nlZ, dnlZ] = ...
        infExact_fast(theta, mean_function, covariance_function, 'likGauss', x, y);
    return;
  end

  % skipping error checks on likelihood, assuming likGauss no matter
  % what the user says

  % determine what needs to be computed
  compute_dalpha = (nargout >= 4);
  compute_dWinv  = (nargout >= 5);
  compute_HnlZ   = (nargout >= 6);

  n = size(x, 1);
  I = eye(n);

  num_cov  = numel(theta.cov);
  num_mean = numel(theta.mean);
  num_hyperparameters = 1 + num_cov + num_mean;

  % initialize output
  dnlZ = theta;

  if (compute_dalpha)
    dalpha.cov  = zeros(n, num_cov);
    dalpha.lik  = zeros(n, 1);
    dalpha.mean = zeros(n, num_mean);
  end

  if (compute_dWinv)
    dWinv.cov  = zeros(n, num_cov);
    dWinv.lik  = zeros(n, 1);
    dWinv.mean = zeros(n, num_mean);
  end

  if (compute_HnlZ)
    HnlZ.value          = zeros(num_hyperparameters);
    HnlZ.covariance_ind = 1:num_cov;
    HnlZ.likelihood_ind = (num_cov + 1);
    HnlZ.mean_ind       = (num_cov + 2):num_hyperparameters;

    % converts to column vector (needed to vectorize A' below)
    vectorize = @(x) (x(:));

    % computes tr(AB)
    product_trace = @(A, B) (vectorize(A')' * B(:));

    % when computing Hessian, we store the following for reuse:
    %
    % - derivatives of mu with respect to mean parameters
    % - V^{-1} K'_i
    dms       = zeros(n, num_mean);
    V_inv_dKs = zeros(n, n, num_cov);
  end

  % convenience handles
  mu = @(varargin) feval(mean_function{:},       theta.mean, ...
                         x, varargin{:});
  K  = @(varargin) feval(covariance_function{:}, theta.cov,  ...
                         x, [], varargin{:});

  % indices of the diagonal entries of an (n x n) matrix
  diag_ind = (1:(n + 1):(n * n))';

  % determine parametrization
  noise_variance = exp(2 * theta.lik);
  high_noise = (noise_variance >= 1e-6);

  % compute posterior if needed
  if (isstruct(y))
    % in case a posterior is provided, we need to calculate:
    %
    % - (y - mu(x))
    % - L = chol(K + sigma^2 I)

    posterior = y;
    alpha = posterior.alpha;

    % derive y from posterior.alpha
    V = K();
    V(diag_ind) = V(diag_ind) + noise_variance;
    y = V * posterior.alpha;

    if (is_chol(posterior.L))
      % high-noise parameterization: posterior.L contains chol(K / sigma^2 + I)

      factor = (1 / noise_variance);
      L = posterior.L;
    else
      % low-noise parameterization: posterior.L contains -inv(K + \sigma^2 I)
      % in this case it's fastest to recompute L

      factor = 1;
      L = chol_robust(V);
    end

    V_inv_times = @(x) solve_chol(L, x) * factor;

  else

    % handle small noise specially to avoid numerical problems (GPML 3.4)
    if (high_noise)
      % high-noise parameterization: posterior.L contains chol(K / sigma^2 + I)

      factor = (1 / noise_variance);

      V = K() * factor;
      V(diag_ind) = V(diag_ind) + 1;

      L = chol_robust(V);
      posterior.L = L;
    else
      % low-noise parameterization: posterior.L contains -(K + \sigma^2 I)^{-1}

      factor = 1;

      V = K();
      V(diag_ind) = V(diag_ind) + noise_variance;

      L = chol_robust(V);
      posterior.L = -solve_chol(L, I);
    end

    V_inv_times = @(x) solve_chol(L, x) * factor;

    y = y - mu();
    alpha = V_inv_times(y);

    posterior.alpha = alpha;
    posterior.sW    = ones(n, 1) * (1 / sqrt(noise_variance));
  end

  % negative log likelihood
  nlZ = sum(log(diag(L))) + ...
        0.5 * (y' * alpha + n * log(2 * pi / factor));

  % calculate (K + \sigma^2 I)^{-1} if needed
  if (is_chol(posterior.L))
    V_inv = V_inv_times(I);
  else
    % desired inverse already calculated in low-noise case
    V_inv = -posterior.L;
  end

  % derivative with respect to log noise scale
  dnlZ.lik = noise_variance * (trace(V_inv) - alpha' * alpha);

  % precompute (K + \sigma^2 I)^{-1}\alpha if necessary; it's used a lot
  if (compute_dalpha || compute_HnlZ)
    V_inv_alpha = V_inv_times(alpha);
  end

  % derivative of alpha with respect to log noise scale
  if (compute_dalpha)
    dalpha.lik = -2 * noise_variance * V_inv_alpha;
  end

  % derivative of diag W^{-1} with respect to log noise scale
  if (compute_dWinv)
    dWinv.lik = 2 * noise_variance * ones(n, 1);
  end

  % second derivative with respect to log noise scale
  if (compute_HnlZ)
    HnlZ.value(HnlZ.likelihood_ind, HnlZ.likelihood_ind) = ...
        2 * noise_variance^2 * ...
        (2 * alpha' * V_inv_alpha - product_trace(V_inv, V_inv)) + ...
        2 * dnlZ.lik;
  end

  % handle gradient/Hessian entries with respect to mean hyperparameters
  for i = 1:num_mean
    dm = mu(i);

    % derivative of nlZ with respect to this mean parameter
    dnlZ.mean(i) = -dm' * alpha;

    % precompute V^{-1} m' if necessary
    if (compute_dalpha || compute_HnlZ)
      V_inv_dm = V_inv_times(dm);
    end

    % derivitive of alpha with respect to this mean parameter
    if (compute_dalpha)
      dalpha.mean(:, i) = -V_inv_dm;
    end

    if (compute_HnlZ)
      dms(:, i) = dm;

      % mean/mean Hessian entries
      for j = 1:i
        d2m_didj = mu(i, j);

        HnlZ.value(HnlZ.mean_ind(i), HnlZ.mean_ind(j)) = ...
            V_inv_dm' * dms(:, j) - ...
            d2m_didj' * alpha;
      end

      % mean/noise Hessian entry
      HnlZ.value(HnlZ.mean_ind(i), HnlZ.likelihood_ind) = ...
          2 * noise_variance * dms(:, i)' * V_inv_alpha;
    end
  end

  % handle gradient/Hessian entries with respect to covariance
  % hyperparameters
  for i = 1:num_cov
    dK = K(i);

    % compute V^{-1} K'; save for Hessian computations if necessary
    V_inv_dK = V_inv_times(dK);
    if (compute_HnlZ)
      V_inv_dKs(:, :, i) = V_inv_dK;
    end

    % derivative of nlZ with respect to this covariance parameter
    dnlZ.cov(i) = 0.5 * (trace(V_inv_dK) - alpha' * dK * alpha);

    % derivative of alpha with respect to this covariance parameter
    if (compute_dalpha)
      dalpha.cov(:, i) = -V_inv_dK * posterior.alpha;
    end

    if (compute_HnlZ)
      % covariance/covariance Hessian entries
      for j = 1:i
        HK = K(i, j);

        HnlZ.value(HnlZ.covariance_ind(i), HnlZ.covariance_ind(j)) = ...
            (y' * V_inv_dKs(:, :, i)) * (V_inv_dKs(:, :, j) * alpha) + ...
            0.5 * (product_trace(V_inv, HK) - ...
                   product_trace(V_inv_dKs(:, :, i), V_inv_dKs(:, :, j)) - ...
                   alpha' * HK * alpha);
      end

      % covariance/mean Hessian entries
      for j = 1:num_mean
        HnlZ.value(HnlZ.mean_ind(j), HnlZ.covariance_ind(i)) = ...
            dms(:, j)' * V_inv_dKs(:, :, i) * alpha;
      end

      % covariance/noise Hessian entry
      HnlZ.value(HnlZ.likelihood_ind, HnlZ.covariance_ind(i)) = ...
          noise_variance * (2 * y' * V_inv_dKs(:, :, i) * V_inv_alpha - ...
              product_trace(V_inv_dKs(:, :, i), V_inv));
    end
  end

  % symmetrize Hessian
  if (compute_HnlZ)
    HnlZ.value = HnlZ.value + tril(HnlZ.value, -1)';
  end

end

function [L,p] = chol_robust(M)
[L,p] = chol(M);                % Try computing Cholesky factor

if p~=0     % Failed Cholesky decomposition, compute nearest SPD matrix
    M = (M + M')/2;         % Ensure M is symmetric
    [U,Sigma,V] = svd(M);
    H = V*Sigma*V';         % Symmetric polar factor H is SPD
    M = (M+H)/2;            
    M = (M + M')/2;         % Ensure symmetry again
    [L,p] = chol(M);        % Retry Cholesky decomposition
    k = 0;
    while p ~= 0            % Failed again, add a small diagonal component
        k = k + 1;
        kappa = 0.05;
        if k > 3; error('Cannot compute Cholesky decomposition.'); end
        mineig = min(eig(M));
        % M = M + (-mineig*k.^2 + eps(mineig))*eye(size(M));
        M = M + kappa*mineig*k.^2*eye(size(M));  % Small nudge
        [L,p] = chol(M);
    end
end

end