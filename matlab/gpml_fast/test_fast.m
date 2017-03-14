D = 10;
N = 100;
x = randn(N, D);
y = randn(N, 1);
TolErr = 1e-12;

mean = {@meanConst}; lik = {@likGauss};
cov = {@covSEard}; cov_fast = {@covSEard_fast};

hyp.mean = 0;
hyp.cov = zeros(eval(feval(cov{:})),1);
hyp.lik = 0.5*log(eps);

%--------------------------------------------------------------------------
display('Testing SE covariance...');

tic
for i = 1:10; [post nlZ dnlZ1] = infExact(hyp, mean, cov, lik, x, y); end
toc

tic
for i = 1:10; [post nlZ dnlZ2] = infExact_fast(hyp, mean, cov_fast, lik, x, y); end
toc

assert(all(abs(dnlZ1.cov - dnlZ2.cov) < TolErr));
rmse = sqrt(sum((dnlZ1.cov(:) - dnlZ2.cov(:)).^2));
display(['Error on each component less than ' num2str(TolErr) '. Total RMSE = ' num2str(rmse) '.']);

%--------------------------------------------------------------------------
display('Testing Matern covariance...');

cov = {@covMaternard,5}; cov_fast = {@covMaternard_fast,5};
hyp.cov = zeros(eval(feval(cov{:})),1);

tic
for i = 1:10; [post nlZ dnlZ1] = infExact(hyp, mean, cov, lik, x, y); end
toc

tic
for i = 1:10; [post nlZ dnlZ2] = infExact_fast(hyp, mean, cov_fast, lik, x, y); end
toc

assert(all(abs(dnlZ1.cov - dnlZ2.cov) < TolErr));
rmse = sqrt(sum((dnlZ1.cov(:) - dnlZ2.cov(:)).^2));
display(['Error on each component less than ' num2str(TolErr) '. Total RMSE = ' num2str(rmse) '.']);

%--------------------------------------------------------------------------
display('Testing periodic covariance...');

cov = {@covPPERard,{@covSEard}}; cov_fast = {@covPPERard_fast,{@covSEard_fast}};
hyp.cov = zeros(eval(feval(cov{:})),1);
hyp.cov(1) = Inf;

tic
for i = 1:10; [post nlZ dnlZ1] = infExact(hyp, mean, cov, lik, x, y); end
toc

tic
for i = 1:10; [post nlZ dnlZ2] = infExact_fast(hyp, mean, cov_fast, lik, x, y); end
toc

assert(all(abs(dnlZ1.cov - dnlZ2.cov) < TolErr));
rmse = sqrt(sum((dnlZ1.cov(:) - dnlZ2.cov(:)).^2));
display(['Error on each component less than ' num2str(TolErr) '. Total RMSE = ' num2str(rmse) '.']);
