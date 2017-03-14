function xs = searchCrossover(x,gpstruct,LB,UB,optimState,options)
%SEARCHCROSSOVER Genetic-algorithm inspired search (crossover operator).

if nargin == 0
    xs = 'xover';
    return;
end

MeshSize = optimState.meshsize;
SearchFactor = optimState.searchfactor;

nvars = length(x);
N = size(optimState.iterList.x,1);
MutationFactor = 1/nvars;

y = bsxfun(@plus, optimState.iterList.fval, ...
    bsxfun(@times,optimState.iterList.fsd,randn(N,options.Nsearch)));

if 0
    [~,index] = min(y,[],1);
    index = index(:);
else
    
    % Fitness scaled by ranking
    for i = 1:options.Nsearch; [~,idx(:,i)] = sort(y(:,i)); end    
    p = 1./((1:N).^0.5);
    p = p ./ sum(p);
    cp = cumsum(p);
    
    r = rand(1,options.Nsearch);
    index = NaN(1,options.Nsearch);
    for i = 1:N
        index(r <= cp(i)) = i;
        r(r <= cp(i)) = Inf;
    end
end

secondParent = optimState.iterList.x(index,:);

% Duplicate first parent
xs = repmat(x, [options.Nsearch,1]);
mutationMatrix = rand(options.Nsearch,nvars) < MutationFactor;
xs(mutationMatrix) = secondParent(mutationMatrix);

end

%--------------------------------------------------------------------------
%MNRND Random sample from the multinomial distribution.
function r = mnrnd(p)
    cdf = cumsum(p,2);
    m = size(p,1);
    samp_k = bsxfun(@gt, cdf, rand(m,1));
    r = diff([zeros(m,1), samp_k],[],2);
end