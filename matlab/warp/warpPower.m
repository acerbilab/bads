% marshalling of parameters and available warping functions
function varargout = warpPower(y,hyp,varargin)
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
          
% Power warping function g(y) and log of the derivative log(g'(y))>0
% or derivatives of the latter w.r.t. ith hyperparameter
function [gy,lgpy] = g(y,hyp,i)
    
  mu = hyp(1);
  lambda = exp(hyp(2));
  
  idx = y > mu;
  
  D = y(idx) - mu + 1;

  if nargin==2                                                 % function values
    gy = y;    
    gy(idx) = (D.^lambda-1)./lambda + mu;
    
    lgpy = zeros(size(y));
    lgpy(idx) = log(D.^(lambda-1));
  else                                                          % derivatives
      
    if i == 1
        gy = zeros(size(y));
        gy(idx) = 1 - D.^(lambda-1);        
        lgpy = zeros(size(y));
        lgpy(idx) = (1 - lambda)./D;
    elseif i == 2
        gy = zeros(size(y));
        gy(idx) = (1 + D.^lambda.*(lambda.*log(D)-1))./lambda;
        lgpy = zeros(size(y));      
        lgpy(idx) = lambda.*log(D);
    end
  end
  
% invert g(y)
function [y,n,d] = ig(z,hyp)
  y = z;
  mu = hyp(1);
  lambda = exp(hyp(2));
  idx = z > mu;
  y(idx) = (1 + lambda.*(z(idx)-mu)).^(1./lambda) + mu - 1;
  n = 0;
  d = 0;