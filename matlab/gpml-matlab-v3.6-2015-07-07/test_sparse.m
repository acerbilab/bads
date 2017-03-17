clear all, close all, s = 25; randn('seed',s), rand('seed',s)
dev = @(x,z) norm(x(:)-z(:))/max([1,norm(x(:)),norm(z(:))]);

n = 194;
x = 5*rand(n,1); y = sin(x); xs = linspace(0,20,1e4)';

mean = {@meanConst}; hyp.mean = 0.2;
sn = 0.2; hyp.lik = log(sn);
ell = 0.2; sf = 2.3; hyp.cov = [log(ell); log(sf)];
v = 0; cov = {@covPPiso,v};

tic
[nlZ,dnlZ,post] = gp(hyp,[],mean,cov,[],x,y);
t = toc;

tic
  K = sparse(feval(cov{:},hyp.cov,x));
  m = feval(mean{:},hyp.mean,x);
  sn2 = sn*sn;
  L = chol(K/sn2+speye(n))';                                % fast column access
  alpha = L'\(L\((y-m)/sn2));
  nlZs = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sn2)/2;
  dnlZs = hyp;

  nz = 2*nnz(L)-n; r = zeros(nz,1); c = zeros(nz,1); v = zeros(nz,1);    % alloc
  k = 0; Z = zeros(n,n);
  for i=n:-1:1
    li = L(:,i); dii = li(i); jj = find(li); li = li(jj); nj = numel(jj);
    z = zeros(nj,1);
    for ii=nj:-1:1
      j = jj(ii);
      if nj>1, zij = li(2:nj)'*Z(jj(2:nj),j); else zij = 0; end
      zij = (i==j)/dii^2 - zij/dii;
      Z(i,j) = zij; Z(j,i) = zij;
      z(ii) = zij;
    end
    r(k+1  ) = i; c(k+1  ) = i;        v(k+1  ) = z(1    ); k = k+1;
    idx = 1:nj-1;
    r(k+idx) = i; c(k+idx) = jj(1+idx); v(k+idx) = z(1+idx); k = k+nj-1;
    c(k+idx) = i; r(k+idx) = jj(1+idx); v(k+idx) = z(1+idx); k = k+nj-1;
  end
  Q = sparse(r,c,v/sn^2,n,n);

%   Z = solve_chol(L',eye(n)); Z = Z.*(K>0);
%   [r,c,v] = find(K); for k=1:numel(r), v(k) = Z(r(k),c(k)); end
%   Q = sparse(r,c,v/sn^2,n,n);

  for i = 1:numel(hyp.cov)
    dK = sparse(feval(cov{:},hyp.cov,x,[],i));
    dnlZs.cov(i) = sum(sum(Q.*dK))/2 - (alpha'*dK*alpha)/2;
  end
  dnlZs.lik = sn2*(trace(Q)-alpha'*alpha);
  for i = 1:numel(hyp.mean)
    dnlZs.mean(i) = -feval(mean{:},hyp.mean,x,i)'*alpha;
  end
ts = toc;

dev(post.alpha,alpha)+dev(nlZ,nlZs) + ...
dev(dnlZ.cov,dnlZs.cov)+dev(dnlZ.lik,dnlZs.lik)+dev(dnlZ.mean,dnlZs.mean)
fprintf('times %1.2f, %1.2f\n',t,ts)