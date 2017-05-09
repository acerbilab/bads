function Bnew = pollBADS2N(B,x,gpstruct,LB,UB,optimState,options)
%POLLBADS2N Poll 2N fixed basis vectors (Generalized Pattern Search).

nvars = numel(x);

if isempty(B)
    Bnew = [eye(nvars); -eye(nvars)];
    Bnew = Bnew*optimState.C';
    
    if 0
        if size(gpstruct.x,1) > (2*nvars + nvars*nvars/4)
            try
                % H = hessian(@(xi_) gppred(xi_,gpstruct), x);
                H = fhess(@(xi_) gppred(xi_,gpstruct), x, [], 'step', optimState.searchmeshsize);
                H = nearestSPD(H,10);
                % eig(H)
                L = chol(H);
                % Bnew = solve_chol(L,Bnew')';
                Hinv = solve_chol(L,eye(nvars));
                Hinv = nearestSPD(Hinv,10);
                [V,D] = eig(Hinv);
                lD = real(log(diag(D)));
                lD(lD < max(lD) + 0.5*log(eps)) = max(lD) + 0.5*log(eps);
                while 1
                    delta = max(lD) - min(lD);
                    if delta < log(1e6); break; end
                    lD = 0.8*lD; 
                end
                D = exp(0.5*lD);
                Bnew = bsxfun(@times, [V'; -V'], [D(:);D(:)]);
                % Hinv = inv(H);
                % Bnew = (chol(Hinv)*Bnew')';
            catch
                % Use diagonal matrix
            end
        end
    end
    
    % Global vector normalization
    D = sqrt(sum(Bnew(1:nvars,:).*Bnew(1:nvars,:),2));
    N = exp(log(D) - mean(log(D)))./D;
    Bnew = bsxfun(@times, Bnew, [N;N]);
    Bnew = bsxfun(@rdivide, Bnew, gpstruct.pollscale);    
    
else
    Bnew = [];
end

end