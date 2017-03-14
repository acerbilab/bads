function es = ESupdate(mu,lambda,iter)

% Create vector of ES weights
es.mu = mu;
es.lambda = lambda;
es.iter = iter;
tot = es.mu + es.lambda;
w = ceil((1./sqrt(1:tot))/sum(1./sqrt(1:tot))*es.lambda);
nonzero = sum(w > 0);
while (sum(w) - es.lambda) > nonzero
    w = max(0, w - 1);
    nonzero = sum(w > 0);
end
delta = sum(w) - es.lambda;
lastnonzero = find(w > 0, 1, 'last');
w(max(1,lastnonzero-delta+1):lastnonzero) = w(max(1,lastnonzero-delta+1):lastnonzero) - 1;

% Create selection mask
cw = cumsum(w) - w + 1;
idx(cw) = 1;
es.selectmask = cumsum(idx(1:end-1));

end