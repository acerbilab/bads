function C = hessianapprox(fun,x0)

nvars = numel(x0);

Cold = eye(nvars);
M = eye(nvars);
n = 2^10;
sigma = 1;
y = zeros(n,1);

% Compute vector weights
mu = 0.5*n;
weights = zeros(1,1,floor(mu));
weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
% weights(1,1,:) = ones(1,mu);
weights = weights./sum(weights);
mueff = 1/sum(weights.^2);
tol = 1e-2;

for i = 1:10
    i
    x = bsxfun(@plus,x0,sigma*randn(n,nvars)*M');
    for j = 1:n; y(j) = fun(x(j,:)); end
    [yord,ord] = sort(y);
    xi = x(ord(1:n/2),:);
    yi = yord(1:n/2);
    
    topx = bsxfun(@minus,xi,x0);            
    C = sum(bsxfun(@times,weights,topx'*topx),3)/sigma^2/n;
    
    if sum((mean(xi,1)-x0).^2) < tol^2
        break
    end
    
    amu = 2;
    cmu = min(1, amu*(mueff-2+1/mueff)/((nvars+2)^2 + amu*mueff/2));
    cmu = 1;
    C = (1-cmu)*Cold + cmu*C;                
    Cold = C;
    
    [V,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
    lambda = diag(D);                        
    lambdared = max(lambda,max(lambda)*1e-12);
    M = real(V)*diag(lambdared);
        
    sigma = sigma/2;    
end




end