D = 10;
N = 100;
x = randn(N, D);
y = randn(N, 1);
TolErr = 1e-12;

mean = {@constant_mean}; lik = {@likGauss};
cov = {@ard_ratquad_covariance}; cov_fast = {@ard_ratquad_covariance_fast};

hyp.mean = 0;
hyp.cov = zeros(eval(feval(cov{:})),1);
hyp.lik = 0.5*log(eps);

%--------------------------------------------------------------------------
display('Testing RQ covariance...');

tic
for i = 1:10; [post,nlZ,dnlZ1,~,~,HnlZ1] = exact_inference(hyp, mean, cov, lik, x, y); end
toc

tic
for i = 1:10; [post,nlZ,dnlZ2,~,~,HnlZ2] = exact_inference_fast(hyp, mean, cov_fast, lik, x, y); end
toc

HnlZ1
HnlZ2

assert(all(abs(dnlZ1.cov - dnlZ2.cov) < TolErr));
rmse = sqrt(sum((dnlZ1.cov(:) - dnlZ2.cov(:)).^2));
Hrmse = sqrt(sum((HnlZ1.value(:) - HnlZ2.value(:)).^2));
display(['Error on each component less than ' num2str(TolErr) '. Total RMSE on gradient = ' num2str(rmse) '. Total RMSE on Hessian = ' num2str(Hrmse)]);