function [hyp,hypw] = gpHyperSVGD(hyp0,gpstruct,options,removeafter)
%GPHYPERSVGD Fit GP hyperparameters via Stein Variational Gradient Descent.
      
if nargin < 3; options = []; end
if nargin < 4 || isempty(removeafter); removeafter = Inf; end

% Perform initial argument checks/transformations
gpstruct = gpset(gpstruct);

nvars = size(gpstruct.x, 2);
Nsamples = options.gpSamples;
%Nsamples = numel(gpstruct.hyp);

hypw = ones(1,Nsamples)/Nsamples;   % Equally weighted samples

% Find fixed parameters (delta or clamped priors), they are not optimized
fixed = []; offset = 0;
for field = {'cov', 'lik', 'mean'}
    num_entries = numel(gpstruct.prior.(field{:}));
    fixed = [fixed,zeros(1,num_entries)];
    for j = 1:num_entries
        if isequal(gpstruct.prior.(field{:}){j}{1},@delta_prior) || ...
                strcmp(gpstruct.prior.(field{:}){j}{1}, 'delta_prior') || ...
                isequal(gpstruct.prior.(field{:}){j}{1},@priorClamped) || ...
                strcmp(gpstruct.prior.(field{:}){j}{1}, 'priorClamped') || ...
                isequal(gpstruct.prior.(field{:}){j}{1},@priorDelta) || ...
                strcmp(gpstruct.prior.(field{:}){j}{1}, 'priorDelta')
            fixed(offset+j) = 1;
        end
    end
    offset = offset + num_entries;      
end
fixed = logical(fixed);

% Get bounds
bounds = unwrap2vec(gpstruct.bounds);
LB = bounds(1:2:end-1); UB = bounds(2:2:end);
lb = LB(~fixed); ub = UB(~fixed);

master_stepsize = 0.1;
h = -1;
auto_corr = 0.9;
method = 'adagrad';
max_iter = options.gpSVGDiters;

dlog_p = @(theta) gp_grad(theta, fixed, gpstruct, Nsamples);

% Rescaling factor (1 for all parameters except GP mean)
rescale = ones(1, size(lb,2));
if ~fixed(end); rescale(end) = std(gpstruct.y); end

% Unwrap starting point, force it to be inside bounds
theta0 = [];
for i = 1:Nsamples
    theta0i = unwrap2vec(hyp0(1+mod(i-1,numel(hyp0))));
    theta0i = theta0i + 0.1.*randn(size(theta0i));
    theta0i = min(max(theta0i(~fixed),lb),ub);
    theta0 = [theta0;theta0i'];
end

try
    % theta0    
    theta = svgd(theta0, dlog_p, max_iter, master_stepsize, h, auto_corr, method, rescale);    
    % theta
    
    kk = numel(theta)/Nsamples;
    for i = 1:Nsamples
        newtheta = unwrap2vec(gpstruct.hyp(1));
        newtheta(~fixed) = min(max(theta(i,:)',lb),ub);
        hyp(i) = rewrap(gpstruct.hyp(1),newtheta);
    end
        
catch
    warning('bads:gpHyperSVGDFail', ['Failed SVGD of hyper-parameters (' num2str(1) ' attempts). GP approximation might be unreliable.']);
    % Return input samples
    for i = 1:Nsamples
        hyp(i) = hyp0(1+mod(i-1,numel(hyp0)));
    end
end


    function dlZ = gp_grad(theta, fixed, gpstruct, Nsamples)
    %GP_OPTIMIZER Wrapper function for GP optimization

        newtheta = unwrap2vec(gpstruct.hyp(1));
        f = @(thetastr) (feval(gpstruct.inf{:}, thetastr, gpstruct.mean, gpstruct.cov, gpstruct.lik, gpstruct.x, gpstruct.y));
        kk = numel(theta)/Nsamples;
        dlZ = [];

        for ii = 1:Nsamples
            %newtheta(~fixed) = theta((1:kk)+(ii-1)*kk);
            newtheta(~fixed) = theta(ii,:);
            thetastruct = rewrap(gpstruct.hyp(1), newtheta);
            [~, ~, dnlZi] = f(thetastruct);
            dnlZi = unwrap2vec(dnlZi);
            dnlZi(fixed) = [];
            dlZ = [dlZ; -dnlZi'];            
        end

    end

end

%--------------------------------------------------------------------------
% The MIT License (MIT)
% 
% The following code is Copyright (c) 2016 Qiang Liu and Dilin Wang
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function  theta = svgd(theta0, dlog_p, max_iter, master_stepsize, h, auto_corr, method, rescale)

%%%%%%%%
% Bayesian Inference via Stein Variational Gradient Descent

% input:
%   -- theta0: initialization of particles, m * d matrix (m is the number of particles, d is the dimension)
%   -- dlog_p: function handle of first order derivative of log p(x)
%   -- max_iter: maximum iterations
%   -- master_stepsize: the general learning rate for adagrad
%   -- h/bandwidth: bandwidth for rbf kernel. Using median trick as default
%   -- auto_corr: momentum term
%   -- method: use adagrad to select the best \epsilon
%   -- rescale: rescaling factor for different dimensions

% output:
%   -- theta: a set of particles that approximates p(x)
%%%%%%%%

if nargin < 4; master_stepsize = 0.1; end;

% for the following parameters, we always use the default settings
if nargin < 5; h = -1; end;
if nargin < 6; auto_corr = 0.9; end;
if nargin < 7; method = 'adagrad'; end;
if nargin < 8; rescale = ones(1,size(theta0,2)); end;

switch lower(method)
    
    case 'adagrad'
        %% AdaGrad with momentum
        theta = theta0;
        
        fudge_factor = 1e-6;
        historial_grad = 0;
        
        for iter = 1:max_iter
            grad = KSD_KL_gradxy(theta, dlog_p, h, rescale);   %\Phi(theta)
            if historial_grad == 0
                historial_grad = historial_grad + grad.^2;
            else
                historial_grad = auto_corr * historial_grad + (1 - auto_corr) * grad.^2;
            end
            adj_grad = grad ./ (fudge_factor + sqrt(historial_grad));
            theta = theta + master_stepsize * adj_grad; % update
        end
        
    otherwise
        error('wrong method');
end
end

function [Akxy, info] = KSD_KL_gradxy(x, dlog_p, h, rescale)
%%%%%%%%%%%%%%%%%%%%%%
% Input:
%    -- x: particles, n*d matrix, where n is the number of particles and d is the dimension of x 
%    -- dlog_p: a function handle, which returns the first order derivative of log p(x), n*d matrix
%    -- h: bandwidth. If h == -1, h is selected by the median trick
%    -- rescale: rescaling factor for different dimensions

% Output:
%    --Akxy: n*d matrix, \Phi(x) is our algorithm, which is a smooth
%    function that characterizes the perturbation direction
%    --info: kernel bandwidth
%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3; h = -1; end % median trick as default

[n, d] = size(x);

%%%%%%%%%%%%%% Main part %%%%%%%%%%
Sqy = dlog_p(x);
Sqy = bsxfun(@times, Sqy, rescale);

% Using rbf kernel as default
x = bsxfun(@rdivide, x, rescale);

XY = x*x';
x2= sum(x.^2, 2);
X2e = repmat(x2, 1, n);

H = (X2e + X2e' - 2*XY); % calculate pairwise distance

% median trick for bandwidth
if h == -1
    h = sqrt(0.5*median(H(:)) / log(n+1));   %rbf_dot has factor two in kernel
end

Kxy = exp(-H/(2*h^2));   % calculate rbf kernel


dxKxy= -Kxy*x;
sumKxy = sum(Kxy,2);
for i = 1:d
    dxKxy(:,i)=dxKxy(:,i) + x(:,i).*sumKxy;
end
dxKxy = dxKxy/h^2;
Akxy = (Kxy*Sqy + dxKxy)/n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.bandwidth = h;

return;
end