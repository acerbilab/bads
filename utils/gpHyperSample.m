function hyps = gpHyperSample(hyp0,gpstruct,options)

covlen = length(gpstruct.hyp(1).cov);
if length(gpstruct.hyp) == 1 || 1
    hyp = gpstruct.hyp(end);
    if NoiseKnown; x0 = hyp.cov;
    else x0 = [hyp.cov;hyp.lik]; end
    if ~isempty(gpstruct.hypmean); x0 = gpstruct.hypmean(:); end

    uncopt = optimset('Display','off','GradObj','on','TolFun',1e-5,'TolX',1e-4,'MaxIter',1);
    samples = fmincon(@(x_) slice_wrapper(x_,gpstruct,-1),x0,[],[],[],[],gpstruct.lb,gpstruct.ub,[],uncopt);

    % hh = hessian(@(x_) slice_wrapper(x_,gpstruct,-1),samples);
    hyp.cov = samples(1:covlen,1);
    if ~NoiseKnown
        hyp.lik = samples(covlen+1:end,1);
        % hyp.lik = max(hyp.lik,-5);
    end
    hyp.w = 1;
else
    hyp = gpstruct.hyp(end);
end

% hyp = gpstruct.hyp;

if NoiseKnown
    x0 = hyp.cov;
else
    x0 = [hyp.cov;hyp.lik];
end
widths = 2*ones(size(x0));

if 1
% vec = combvec([-1 0 1],[-1 0 1],[-1 0 1])';
    x0 = samples;
    vec = combvec([-2 -1 0 1 2],[-2 -1 0 1 2],[-2 -1 0 1 2])';
    x = bsxfun(@plus,x0',vec);

    if 0
        for d = 1:length(x0)
            xx = linspace(x0(d)-widths(d),x0(d)+widths(d),21);
            for i = 1:length(xx)
                xt = x0; xt(d) = xx(i);
                lZ(i) = slice_wrapper(xt,gpstruct);
            end
            subplot(1,length(x0),d);
            plot(xx,lZ);
            ylim([max(lZ)-10,max(lZ)]);
            drawnow;
        end
    end

    for i = 1:size(x,1)
        lZ(i) = slice_wrapper(x(i,:),gpstruct);                    
    end
    lZ(isnan(lZ)) = -Inf;
    nz = max(lZ);
    hypw = exp(lZ-nz)./sum(exp(lZ-nz));
    gpstruct.hypmean = sum(bsxfun(@times,hypw(:),x),1);
    %hypw
    %x
    %gpstruct.hypmean

    [hypword,index] = sort(hypw,'descend');
    samples = x(index(1:Nsamples),:)';


    hypword = hypword(1:Nsamples)/sum(hypword(1:Nsamples));                
    cumw = 1 - cumsum(hypword);
    lastgood = find(cumw < 0.01,1);
    samples
    cumw

    if ~isempty(lastgood); cumw(lastgood+1:end) = []; end

    Nsamples = length(cumw);
    hypword(Nsamples+1:end) = [];
    samples(:,Nsamples+1:end) = [];
    hypword = hypword/sum(hypword);
    gpstruct.hyp(Nsamples+1:end) = [];

else
    word = [];
    while 1
        if 1
            uncopt = optimset('Display','off','GradObj','on','TolFun',0.01,'TolX',1e-5,'MaxFunEval',200);
            t0 = fmincon(@(x_) slice_wrapper(x_,gpstruct,-1),x0,[],[],[],[],gpstruct.lb,gpstruct.ub,[],uncopt);

            H = hessian(@(x_) slice_wrapper(x_,gpstruct,-1), t0);
            C = inv(H);
            try
                samples = mvnrnd(ones(Nsamples,1)*t0',C)';

                for i = 1:Nsamples
                    hypword(i) = slice_wrapper(samples(:,i),gpstruct,-1)/mvnpdf(samples(:,i),t0,C);
                end                        
                hypword = hypword./sum(hypword);

            catch
                C
                samples = t0(:);
                % samples = mvnrnd(ones(Nsamples,1)*t0',C.*eye(size(C,1)));
                hypword = 1;
            end


            hypword

        %try
            %ncalls = 0;
            %samples = slice_sample(Nsamples,Nsamples*5,5,@slice_wrapper,x0,widths,0,0,gpstruct);
            break;
        %catch
        %    warning('Error in sampling... trying to recover.');
        %    x0 = zeros(1,length(x0));
        %end
        else

            gpstruct.prior.cov = [];
            gpstruct.prior.lik = [];

            uncopt = optimset('Display','off','GradObj','on','TolFun',0.1,'TolX',1e-4,'MaxFunEval',100);
            t0 = fmincon(@(x_) slice_wrapper(x_,gpstruct,-1),x0,[],[],[],[],gpstruct.lb,gpstruct.ub,[],uncopt);

            chol_Sigma = 4*eye(length(x0));
            samples = zeros(length(x0),1+Nsamples*100);
            samples(:,1) = x0;
            cur_log_like = [];
            for i = 2:size(samples,2)
                [samples(:,i), cur_log_like] = gppu_elliptical(samples(:,i-1), chol_Sigma, @(x_) slice_wrapper(x_,gpstruct,1), cur_log_like);
            end
            samples = samples(:,end-Nsamples*20+1:20:end);
            break;
            % error('a');


        end
    end
end
% samples
% ncalls

    function [lZ,dlZ] = slice_wrapper(x,gpstruct,sig)
        ncalls = ncalls + 1;
        if nargin < 3; sig = 1; end
        hyp.cov = x(1:covlen);
        if ~NoiseKnown; hyp.lik = x(covlen+1); end
        if ~fixedMean; hyp.mean = x(end); end
        % hyp.lik = max(hyp.lik,-5);
        
        if nargout == 1
            nlZ = gp(hyp,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);            
        else
            [nlZ,ndlZstruct] = gp(hyp,gpstruct.inf,gpstruct.meanfunc,gpstruct.covfunc,gpstruct.likfunc,gpstruct.x,gpstruct.y);
            
            dlZ = ndlZstruct.cov(:);
            if ~NoiseKnown; dlZ = [dlZ; ndlZstruct.lik(:)]; end
            if ~fixedMean; dlZ = [dlZ; ndlZstruct.mean(:)]; end
            dlZ = sig*(-dlZ);
        end
        lZ = sig*(-nlZ);
        
    end


end