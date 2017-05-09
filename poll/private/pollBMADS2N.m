function Bnew = pollBMADS2N(B,x,gpstruct,LB,UB,optimState,options)
%POLLBMADS2N Poll 2N fixed basis vectors (Generalized Pattern Search).

nvars = numel(x);

if isempty(B)
        % Axis-oriented poll vector set
        u = [eye(nvars); -eye(nvars)];
        B = eye(nvars);
        
        if optimState.iter > 1 || 1
            [q,~] = qr(randn(nvars));
            % u = randn(nvars,nvars);
            % u = bsxfun(@rdivide,u,sqrt(sum(u.*u,2)));
            Bnew = q*B;
        elseif size(gpstruct.x,1) > 2*nvars && 0

            % Add small jitter to prevent singularity
            xscatter = gpstruct.x + 0.1*bsxfun(@times,MeshSize,randn(size(gpstruct.x)));
            xscatter = bsxfun(@times,xscatter,optimState.scale);

            % Compute vector weights
            mu = size(gpstruct.x,1)/2;
            weights = zeros(1,1,floor(mu));
            weights(1,1,:) = log(mu+1/2)-log(1:floor(mu));
            weights = weights./sum(weights);

            % Compute top vectors
            [~,index] = sort(gpstruct.y,'ascend');

            % Compute weighted covariance matrix wrt X0
            topx = bsxfun(@minus,xscatter(index(1:floor(mu)),:),x);
            C = sum(bsxfun(@times,weights,topx'*topx),3);

            [PositiveBasis,e] = eig(C);
            u = [PositiveBasis; -PositiveBasis];
        else
            if mod(optimState.pollcount,2) == 0 && 1
                % Axis-oriented poll vector set
                u = [eye(nvars); -eye(nvars)];
                PositiveBasis = eye(nvars);
            elseif 0
                % Mixed axis poll vector set
                PositiveBasis = zeros(nvars,nvars);
                idx = 1;
                for i = 1:nvars
                    while 1
                        j = randi(nvars-1); if j == i; j = j + 1; end
                        CandidateVector = zeros(1, nvars);
                        CandidateVector(i) = 1;
                        CandidateVector(j) = 0.25*2*(randi(2)-1.5);
                        if ~any(all(bsxfun(@eq, CandidateVector, PositiveBasis),2)) ...
                                && ~any(all(bsxfun(@eq, -CandidateVector, PositiveBasis),2))
                            PositiveBasis(i,:) = CandidateVector;
                        end
                        break;
                    end
                end
                PositiveBasis = PositiveBasis;
                u = [PositiveBasis; -PositiveBasis];
            end
        end
        
elseif 0
        
        v = ((xnew - x)./optimState.scale)./MeshSize;

        if ~isempty(PollBasis)
            t1 = sum(abs(bsxfun(@minus,v,PollBasis)),2); 
            t2 = sum(abs(bsxfun(@minus,-v,PollBasis)),2);
            if all(t1 > optimState.TolMesh/2) && all(t2 > optimState.TolMesh/2)
                PositiveBasis = [PositiveBasis; v];
            end
        end

        % Store polled vectors in this iteration
        PollBasis = [PollBasis; v];

        % Generate new MADS set orthogonal to existing basis
        u = randn(options.Nbasis,nvars);
        u = bsxfun(@rdivide,u,sqrt(sum(u.*u,2)));

        for j = 1:size(PositiveBasis,1)
            u = u - (u*PositiveBasis(j,:)'./(PositiveBasis(j,:)*PositiveBasis(j,:)'))*PositiveBasis(j,:);
        end
        u = bsxfun(@rdivide,u,sqrt(sum(u.*u,2)));
        u = [u; -PositiveBasis];

        xpoll = bsxfun(@plus,x,bsxfun(@times,u*MeshSize,optimState.scale));

        % Remove poll vectors outside bounds
        rf = logical(sum(bsxfun(@gt,xpoll,UB) | bsxfun(@lt,xpoll,LB),2));
        xpoll(rf,:) = [];

        % Remove previously evaluted vectors
        % xpoll = setdiffrows(xpoll,funlog.X);
else
    Bnew = [];
                    
end


end