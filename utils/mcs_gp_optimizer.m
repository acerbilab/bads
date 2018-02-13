function [nlZ, dnlZ, HnlZ] = mcs_gp_optimizer(gpstruct, theta)
%MCS_GP_OPTIMIZER Wrapper function for GP optimization

fixed = gpstruct.fixed;
hyp = gpstruct.hyp;
inf = gpstruct.inf; 
mean = gpstruct.mean; 
cov = gpstruct.cov; 
lik = gpstruct.lik; 
x = gpstruct.x;
y = gpstruct.y;

newtheta = unwrap2vec(hyp);
newtheta(~fixed) = theta;
theta = rewrap(hyp, newtheta);

f = @() (feval(inf{:}, theta, mean, cov, lik, x, y));

if nargout <= 1
    [~, nlZ] = f();
    return;

elseif nargout == 2
    [~, nlZ, dnlZ] = f();

elseif nargout > 2
    [~, nlZ, dnlZ, ~, ~, HnlZ] = f();

    HnlZ = HnlZ.value;
    HnlZ(fixed,:) = [];
    HnlZ(:,fixed) = [];
end

dnlZ = unwrap2vec(dnlZ);
dnlZ(fixed) = [];

end