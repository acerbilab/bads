function [K,Mx,xe] = covGrid(cov, xg, hyp, x, z, i)

% covGrid - Kronecker covariance function based on a grid.
%
% The grid g is represented by its p axes xg = {x1,x2,..xp}. An axis xi is of
% size (ni,di) and the grid g has size (n1,n2,..,np,D), where D=d1+d2+..+dp.
% Hence, the grid contains N=n1*n2*..*np data points. The axes do neither need
% to be sorted in any way nor do they need to be 1d i.e. ni>=1.
%
% The covGrid function can be used to expand the cell array xg into an expanded
% multivariate grid xe of size (N,D) via:
%     [xe,nx,Dx]  = covGrid('expand',xg);                             => mode 1)
% The operation can be reverted by:
%     xg = covGrid('factor',{xe,ng,Dg});                              => mode 2)
%
% Given scattered data x of size (n,D), we can create a grid xg covering
% the support of x using:
%     xg = covGrid('create',x,equi,k);                                => mode 3)
% The flag equi (default value 1) can be used to enforce an equispaced
% grid. The integer k, indicates the number of grid points per dimension. If
% k is a real number from (0,1], then the number of grid points equals
% k*numel(unique(x(:,1))).
% We require at least two different components per dimension.
%
% The variables v={x,z} can either be a) grid indices or b) data points.
% a) The variable v has size (nv,1) and contains integers from [1,N]. Then
%    the datapoints are obtained as g2 = reshape(g,N,D); v = g2(v,:).
% b) The variable v has size (nv,D) and directly represents the data points.
% At the moment, this is only supported for z.
%
% An arbitrary data point x -- be it an index vector of size (n,1) or a data
% point of size (n,D) -- is converted into a regular data point xx of 
% size (n,D) by:
%     xx = covGrid('idx2dat',xg,x);                                   => mode 4)
% If x is already of size (n,D), xx will simply equal x.
% 
% The resulting covariance matrix is given by:
%   K = kron( kron(...,K{2}), K{1} ) = K_p x .. x K_2 x K_1.
%
% The hyperparameters are:
% hyp = [ hyp_1
%         hyp_2
%          ..
%         hyp_p ],
%
% Copyright (c) by Hannes Nickisch and Andrew Wilson 2015-03-05.
%
% See also COVFUNCTIONS.M, INFGRID.M.

if nargin<2, error('Not enough parameters provided.'), end

% mode 1) expand axes xg representation into full grid x
if     strcmp(cov,'expand')          % call: [xe,nx,Dx]  = covGrid('expand',xg);
  [K,Mx,xe] = expandgrid(xg); return

% mode 2) factor full x grid into axes representation xg
elseif strcmp(cov,'factor')           % call: xg = covGrid('factor',{xe,ng,Dg});
  K = factorgrid(xg{:}); return

% mode 3) create axes representation xg from scattered data
elseif strcmp(cov,'create')             % call: xg = covGrid('create',x,equi,k);
  if nargin<3, equi = 1; else equi = hyp; end               % set default values
  if nargin<4, k = 1; else k = x; end, x = xg;                % set input params
  p = size(x,2); xg = cell(p,1);                               % allocate result
  if numel(k)>0, k = ones(p,1).*k(:); end              % enforce vector-valued k
  for j=1:p                                            % iterate over dimensions
    u = sort(unique(x(:,j))); if numel(u)<2, error('Two few unique points.'),end
    if isempty(k)                              % determine number of grid points
      if equi
        ngj = ceil( (u(end)-u(1))/min(abs(diff(u))) );     % use minimum spacing
      else
        ngj = numel(u);
      end
    elseif 0<=k(j) && k(j)<=1
      ngj = ceil(k(j)*numel(u));
    else
      ngj = k(j);
    end
    du = (u(end)-u(1))/ngj; bu = [u(1)-5*du, u(end)+5*du];
    if equi                                                    % equispaced grid
      xg{j} = linspace(bu(1),bu(2),max(ngj,5))';        % at least 5 grid points
    else                                                   % non-equispaced grid
      [idx,xgj] = kmeans(u,min(numel(u),ngj-2)); xgj = sort(xgj(:))';  % cluster
      nb = ngj-numel(xgj); nb1 = floor(nb/2); nb2 = nb - nb1; % size of boundary
      xg1 = linspace(bu(1),xgj(1),nb1+1); xg2 = linspace(xgj(end),bu(2),nb2+1);
      xg{j} = [xg1(1:nb1), xgj, xg2(1+(1:nb2))]';
    end
  end
  K = xg; return

% mode 4) convert possible index vector into data space
elseif strcmp(cov,'idx2dat')               % call: xx = covGrid('idx2dat',xg,x);
  N = 1; for i=1:numel(xg), N = N*size(xg{i},1); end
  if isidx(hyp,N), xe = expandgrid(xg); K = xe(hyp,:); else K = hyp; end, return
end

% mode 0) regular covariance function computations
p = numel(xg); ng = zeros(p,1); Dg = zeros(p,1);   % number of Kronecker factors
if numel(cov)~=p, error('We require p factors.'), end
for ii = 1:p                                 % iterate over covariance functions
  [ng(ii),Dg(ii)] = size(xg{ii});
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  D = Dg(ii); j(ii) = cellstr(num2str(eval(feval(f{:}))));  % collect nbr hypers
end

if nargin<4                                        % report number of parameters
  K = char(j(1)); for ii=2:length(cov), K = [K, '+', char(j(ii))]; end, return
end
if nargin<5, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

v = [];               % v vector indicates to which covariance parameters belong
for ii = 1:p, v = [v repmat(ii, 1, eval(char(j(ii))))]; end
if nargin==6 && i>length(v), error('Unknown hyperparameter'), end

equi = true(p,1);     % grid along dimension is equispaced -> Toeplitz structure
for ii=1:p
  if ng(ii)>1       % diagnose Toeplitz structure if data is linearly increasing
    dev = abs(diff(xg{ii})-ones(ng(ii)-1,1)*(xg{ii}(2,:)-xg{ii}(1,:)));
    equi(ii) = max(dev(:))<1e-9;
  end
end

N = prod(ng); n = size(x,1); D = sum(Dg);     % expanded grid and data dimension
ix = isidx(x,N);               % determine whether x is an index or a data array
if ~ix && nargout>1                                     % off-grid interpolation
  error('No off-grid values for x allowed.')
end
if dg               % evaluate as full dense vector for diagonal covariance case
  K = 1;                       % xg is not assumed to form a grid for z = 'diag'
  for ii = 1:length(cov)                       % iteration over factor functions
    f = cov(ii); if iscell(f{:}), f = f{:}; end % expand cell array if necessary
    d = sum(Dg(1:ii-1))+(1:Dg(ii));                     % dimensions of interest
    if nargin<6, i = 0; vi = 0; else vi = v(i); end; % which covariance function
    if i<=length(v)
      if ix, xii = xg{ii}; else xii = x(:,d); end  % switch Kronecker/plain prod
      if ii==vi
        j = sum(v(1:i)==vi);                % which parameter in that covariance
        Kj = feval(f{:}, hyp(v==ii), xii, z, j);        % deriv Kronecker factor
      else
        Kj = feval(f{:}, hyp(v==ii), xii, z);           % plain Kronecker factor
      end
      if ix, K = kron(K,Kj); else K = K.*Kj; end   % switch Kronecker/plain prod
    else error('Unknown hyperparameter')
    end
  end
  if ix, K = K(x); end, return
end

if isidx(z,N), iz = z; z = covGrid('expand',xg); z = z(iz,:); end     % expand z
K = cell(p,1);                                    % covariance Kronecker factors
for ii = 1:length(cov)                         % iteration over factor functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  d = sum(Dg(1:ii-1))+(1:Dg(ii));                       % dimensions of interest
  if isnumeric(z) && ~isempty(z)                                   % cross terms
    zd = z(:,d);
  elseif xeqz && equi(ii)                           % we have Toeplitz structure
    zd = xg{ii}(1,:);
  else                                                        % symmetric matrix
    zd = z;
  end
  if nargin<6, i = 0; vi = 0; else vi = v(i); end;   % which covariance function
  if i<=length(v)
    if ii==vi
      j = sum(v(1:i)==vi);                  % which parameter in that covariance
      K{ii} = feval(f{:}, hyp(v==ii), xg{ii}, zd, j);   % deriv Kronecker factor
    else
      K{ii} = feval(f{:}, hyp(v==ii), xg{ii}, zd);      % plain Kronecker factor
    end
  else error('Unknown hyperparameter')
  end
  if xeqz && equi(ii), K{ii} = {'toep',K{ii}}; end      % make Toeplitz explicit
end

if ~xeqz                                                    % expand cross terms
  Ks = K; K = Ks{1}; for ii = 2:p, K = kron1(Ks{ii},K); end
  if ix, if numel(x)~=N || max(abs(x-(1:N)'))>0, K = K(x,:); end
  else   K = Mx*K;
  end
end
if nargout>1, if ix, Mx = sparse(1:n,x,1,n,N); end, end
if nargout>2, xe = covGrid('expand',xg); end

% perform kron along first dimension only
% the code is equivalent to the following loop
%   z = zeros(size(x,1)*size(y,1),size(x,2));
%   for i=1:size(z,2), z(:,i) = kron(x(:,i),y(:,i)); end
function z = kron1(x,y)
  nx = size(x,1); ny = size(y,1);
  z = repmat(reshape(x,1,nx,[]),[ny,1,1]).*repmat(reshape(y,ny,1,[]),[1,nx,1]);
  z = reshape(z,nx*ny,[]);

function r = isidx(i,N)     % check whether i represents an integer index vector
  r = false;
  if numel(i)>0 && ~strcmp(i,'diag') && size(i,2)==1 && ndims(i)==2
    if max(abs(i-floor(i)))<1e-13
      if 0<min(i) && max(i)<=N, r = true; end
    end
  end

function [x,ng,Dg] = expandgrid(xg)                    % expand a Kronecker grid
  p = numel(xg); x = xg{1};                                 % expanded grid data
  ng = zeros(p,1); Dg = zeros(p,1); [ng(1),Dg(1)] = size(xg{1});
  for i=2:p
    szx = size(x); [ng(i),Dg(i)] = size(xg{i});
    xold = repmat(reshape(x,[],1,szx(end)),[1,ng(i),1]);
    xnew = repmat(reshape(xg{i},[1,ng(i),Dg(i)]),[size(xold,1),1,1]);
    x = reshape(cat(3,xold,xnew),[szx(1:end-1),ng(i),szx(end)+Dg(i)]);
  end
  x = reshape(x,[],size(x,ndims(x)));

function xg = factorgrid(x,ng,Dg)                      % factor a Kronecker grid
 p = numel(ng); xg = cell(p,1);          % extract individual grid components xg
for i=1:p
  x = reshape(x,[prod(ng(1:i-1)), ng(i), prod(ng(i+1:end)), sum(Dg)]);
  xg{i} = reshape(x(1,:,1,sum(Dg(1:i-1))+(1:Dg(i))), ng(i), Dg(i));
end