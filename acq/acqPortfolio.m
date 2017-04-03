function [hedge, acqIndex, ymu, ys, fmu, fs] = acqPortfolio(method, hedge, u, f, fs, gpstruct, optimState, options, SufficientImprovement, fvalold, MeshSize)
%ACQPORTFOLIO Evaluate and update portfolio of acquisition functions (unsupported).

if nargin < 9; fvalold = []; end
if nargin < 10; MeshSize = []; end

switch lower(method(1:3))
    case 'acq'
        
        [acqIndex,~,ymu,ys,fmu,fs] = ...
            acqHedge(u,f,fs,gpstruct,optimState,options,SufficientImprovement);

        % Initialize hedge struct
        if isempty(hedge)
            [hedge.g,hedge.str] = acqHedge();
            hedge.n = numel(hedge.g);
            hedge.count = 0;
            hedge.gamma = options.HedgeGamma;
            hedge.beta = options.HedgeBeta;
            hedge.decay = options.HedgeDecay;
        end

        hedge.count = hedge.count + 1;
        % gammaHedge = min(1/nHedge, sqrt(log(nHedge)/(nHedge*hedge.count)));

        hedge.p = exp(hedge.beta*(hedge.g - max(hedge.g)))./sum(exp(hedge.beta*(hedge.g - max(hedge.g))));
        hedge.p = hedge.p*(1-hedge.n*hedge.gamma) + hedge.gamma;

        hedge.chosen = find(rand() < cumsum(hedge.p),1);
        if hedge.gamma == 0
            hedge.phat = ones(size(hedge.g));
        else
            hedge.phat = Inf(size(hedge.p));
            hedge.phat(hedge.chosen) = hedge.p(hedge.chosen);
        end
        
    case 'upd'

        for iHedge = 1:hedge.n
            uHedge = u(min(iHedge,end),:);

            %gpstructnew = gpTrainingSet(gpstructnew, ...
            %    options.gpMethod, ...
            %    uHedge, ...
            %    [], ...
            %    optimState, ...
            %    options, ...
            %    0);

            if iHedge == hedge.chosen
                fHedge = f; fsHedge = fs;
            elseif hedge.gamma == 0
                % Compute estimated function value at point
                [~,~,fHedge,fs2Hedge] = gppred(uHedge,gpstruct);
                if numel(gpstructnew.hyp) > 1
                    fHedge = weightedsum(gpstruct.hypweight,fHedge,1);
                    fs2Hedge = weightedsum(gpstruct.hypweight,fs2Hedge,1);
                end
                fsHedge = sqrt(fs2Hedge);
            else
                fHedge = 0; fsHedge = 1;
            end
            
            if fsHedge == 0
                er = max(0, fvalold - fHedge);
            elseif isfinite(fHedge) && isfinite(fsHedge) && isreal(fsHedge) && fsHedge > 0
                % Probability of improvement
                gammaz = (fvalold - fHedge)./fsHedge;
                fpi = 0.5*erfc(-gammaz/sqrt(2));            

                % Expected reward
                er = fsHedge.*(gammaz.*fpi + exp(-0.5*(gammaz.^2))/sqrt(2*pi));
            else
                er = 0;
            end

            hedge.g(iHedge) = hedge.decay*hedge.g(iHedge) + er/hedge.phat(iHedge)/MeshSize;
        end
         
end