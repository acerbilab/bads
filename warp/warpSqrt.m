% marshalling of parameters and available warping functions
function varargout = warpSqrt(y,hyp,varargin)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  
  % return number of parameters
  m = 2;
  if nargin<1 && nargout>0, varargout{1} = m; return, end
  
  if nargin <= 2
      [varargout{:}] = g(y,hyp);
  elseif nargin == 3
      if strcmpi(varargin{1},'inv')
          [varargout{:}] = ig(y,hyp);
      elseif isnumeric(varargin{1})
          [varargout{:}] = g(y,hyp,varargin{1});
      end
  else
      error('Derivative of inverse not supported yet.');
  end
          
% Square root warping function g(y) and log of the derivative log(g'(y))>0
% or derivatives of the latter w.r.t. ith hyperparameter
function [gy,lgpy] = g(y,hyp,i)
    
  mu = hyp(1);
  delta = exp(hyp(2));
  
  idx = y > mu;

  if nargin==2                                                 % function values
    gy = y;    
    gy(idx) = delta*sqrt((y(idx) - mu)/delta +1) + mu - delta;
    
    lgpy = zeros(size(y));
    lgpy(idx) = 0.5*log(delta) - 0.5*log(y(idx)-mu+delta) - log(2);
  else                                                          % derivatives
      
    if i == 1
        gy = zeros(size(y));
        gy(idx) = 1 - 0.5./sqrt(1 + (y(idx) - mu)./delta);        
        lgpy = zeros(size(y));
        lgpy(idx) = 0.5./(y(idx)-mu+delta);      
    elseif i == 2
        gy = zeros(size(y));
        gy(idx) = 0.5*(y(idx) - 2*delta.*(-1 + sqrt(1 + (y(idx) - mu)./delta)) - mu)./sqrt(1 + (y(idx) - mu)./delta);
        lgpy = zeros(size(y));      
        lgpy(idx) = 0.5*(y(idx)-mu)./(y(idx)-mu+delta);
    end
  end
  
% invert g(y)
function [y,n,d] = ig(z,hyp)
  y = z;
  mu = hyp(1);
  delta = exp(hyp(2));
  idx = z > mu;
  y(idx) = delta*(((z(idx)-mu +delta)./delta).^2-1) + mu;
  n = 0;
  d = 0;