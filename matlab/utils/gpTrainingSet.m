function gpstruct = gpTrainingSet(gpstruct,method,uc,ui,optimState,options,retrain)
%GPTRAININGSET Define training input set of Gaussian Process

logindex = 1:optimState.Xmax;
U = optimState.U(logindex,:);
Y = optimState.Y(logindex);
if isfield(optimState,'S'); S = optimState.S(logindex); end
D = numel(uc);
rotate_gp = isfield(gpstruct,'C') && ~isempty(gpstruct.C);            

switch(lower(method))
    case 'add'
        if options.FitnessShaping
            ystar = gpstruct.nonlinf(ui(1), gpstruct.nonlinmu, gpstruct.deltay);
        else
            ystar = ui(1);
        end
        try
            if options.NoiseObj; ysd = ui(2); else ysd = []; end
            gpstruct.post = update_posterior(gpstruct.hyp(1), gpstruct.mean, ...
                  gpstruct.cov, gpstruct.x, gpstruct.post, uc, ystar);          
            if rotate_gp; uc = uc*gpstruct.C'; end     
            gpstruct.x = [gpstruct.x; uc];
            gpstruct.y = [gpstruct.y; ystar];
            if options.NoiseObj; gpstruct.sd = [gpstruct.sd; ysd]; end
        catch
            % Posterior update failed, point was not added
        end
        return;
          
    case 'local'
        
        if isempty(ui)
            % Distance from the central point
            if all(gpstruct.lenscale == 1) || ~options.ScaleDistance % Not supported
                dist = sum(bsxfun(@minus,U,uc).^2,2);
            else
                dist = sum(bsxfun(@rdivide,bsxfun(@minus,U,uc),gpstruct.lenscale).^2,2);
            end
        else
            % Distance between vector and set of poll vectors
            dist = zeros(size(U,1),size(ui,1));
            if all(gpstruct.lenscale == 1) || ~options.ScaleDistance
                for i = 1:size(ui,1)
                    dist(:,i) = sum(bsxfun(@minus,U,ui(i,:)).^2,2);
                end
            else
                for i = 1:size(ui,1)
                    dist(:,i) = sum(bsxfun(@rdivide,bsxfun(@minus,U,ui(i,:)),gpstruct.lenscale).^2,2);
                end
            end
            dist = min(dist,[],2);
        end
        [distord,ord] = sort(dist,'ascend');
        ntrain = min(optimState.Xmax,options.Ndata);
        % ntrain = min(ntrain, sum(distord <= 4))
        
        if 1
            index = 1:ntrain;
        else
            nclose = find([0;distord/optimState.MeshSize^2] < 1,1,'last') - 1;
            if nclose > ntrain
                y = Y(ord(1:nclose));
                [y,ordy] = sort(y,'ascend');
                index = ordy(1:options.Ndata);
            else
                index = 1:ntrain;
            end
        end
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index));
        
    case 'cluster'
                
        if isempty(ui)
            % Distance from the central point
            dist = udist(U,uc,gpstruct.lenscale,optimState);
        else
            % Distance between vector and set of poll vectors
            dist = zeros(size(U,1),size(ui,1));
            for i = 1:size(ui,1)
                dist(:,i) = udist(U,ui(i,:),gpstruct.lenscale,optimState);
            end
            dist = min(dist,[],2);
        end
        [distord,ord] = sort(dist,'ascend');
        
        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*gpstruct.effectiveradius;
        ntrain = max(options.MinNdata,max(options.Ndata-options.BufferNdata, min(options.Ndata, sum(distord <= radius^2))));
        ntrain = min(ntrain, optimState.Xmax);
        % ntrain = max(floor(ntrain/2),min(ntrain, sum(distord <= radius^2)));
        
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Cluster observations
        index = 1:ntrain;
        
        if ntrain > options.Ndata && 0
            % Clustering does not account for periodic data (unused though)
            index = 1:ceil(min(optimState.Xmax,options.Ndata)/2);
            gpstruct.x = U(ord(index),:);
            gpstruct.y = Y(ord(index),:);        
            
            index = size(gpstruct.x,1)+1:ntrain;            
            x = U(ord(index),:);
            y = Y(ord(index),:);
            idx = kmeans([x,y],floor(options.Ndata/2));
            nClusters = max(idx);
            gpstruct.x = zeros(nClusters,D);
            gpstruct.y = zeros(nClusters,1);
            for iCluster = 1:nClusters
                subset = idx == iCluster;
                [ymin,pos] = min(y(subset));
                gpstruct.y(end+1) = ymin;
                xx = x(subset,:);
                gpstruct.x(end+1,:) = xx(pos,:); 
            end
        end
        
        nvars = size(U,2);
        % Add safeguarded points
        for d = 1:nvars
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);        
        
    case 'clusternew'
                
        if isempty(ui)
            % Distance from the central point
            dist = udist(U,uc,gpstruct.lenscale,optimState);
        else
            % Distance between vector and set of poll vectors
            dist = zeros(size(U,1),size(ui,1));
            for i = 1:size(ui,1)
                dist(:,i) = udist(U,ui(i,:),gpstruct.lenscale,optimState);
            end
            dist = min(dist,[],2);
        end
        [distord,ord] = sort(dist,'ascend');
        
        % The first set only contains nearest-neighbors (up to any distance)
        n1 = min([options.MinNdata, optimState.Xmax]);
        
        % The second set is filled with remaining data up to radius
        radius = options.gpRadius*gpstruct.effectiveradius;
        nclose = sum(distord <= radius^2);
        n2 = min([options.Ndata, nclose, optimState.Xmax]);
                
        % Cluster observations        
        if options.gpCluster && n1 < n2 && n2 < nclose
            % Clustering does not account for periodic data (unused though)
            index = 1:n1;            
            nc = n2 - n1;
            nfree = nclose - n1;
            
            clindex = n1+1:nclose;            
            x = U(ord(clindex),:);
            y = Y(ord(clindex),:);
            idx = kmeans([x,y],nc);
            nClusters = max(idx);
            for iCluster = 1:nClusters
                subset = idx == iCluster;
                pos = randi(numel(subset));
                % [~,pos] = min(y(subset));
                index = [index,n1+subset(pos)];
            end
        else
            index = 1:max(n1,n2);
        end
        
        nvars = size(U,2);
        % Add safeguarded points
        for d = 1:nvars
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
        index = unique(index);
        % nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);
        
        numel(index)
        
        
    case 'clustersafe'
                
        if isempty(ui)
            % Distance from the central point
            dist = udist(U,uc,gpstruct.lenscale,optimState);
        else
            % Distance between vector and set of poll vectors
            dist = zeros(size(U,1),size(ui,1));
            for i = 1:size(ui,1)
                dist(:,i) = udist(U,ui(i,:),gpstruct.lenscale,optimState);
            end
            dist = min(dist,[],2);
        end
        [distord,ord] = sort(dist,'ascend');
        
        if isempty(ui)
            % Distance from the central point
            dist2 = udist(U(ord,:),uc,1,optimState);
        else
            % Distance between vector and set of poll vectors
            dist2 = zeros(size(U,1),size(ui,1));
            for i = 1:size(ui,1)
                dist2(:,i) = udist(U(ord,:),ui(i,:),1,optimState);
            end
            dist2 = min(dist2,[],2);
        end
        [distord2,ord2] = sort(dist2,'ascend');
        
        % Add closest points according to gp length scale
        index = 1:min(options.MinNdata,optimState.Xmax);
        
        % Add closest points within double mesh scale unit
        
        % Keep only points within a certain (rescaled) radius from target
        radius = 2;
        nclose = sum(distord2 <= radius^2);
        extra = setdiff(ord2(1:nclose), index)';
        index = [index, extra(1:min(end,options.BufferNdata))];
        ntrain = numel(index);
                        
        nvars = size(U,2);
        % Add safeguarded points
        for d = 1:nvars
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
        index = unique(index);
        nextra = numel(index) - ntrain;
        if nextra > 0; nextra, end
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);        
        
    case 'clusterextra'
                
        if isempty(ui)
            % Distance from the central point
            dist = udist(U,uc,gpstruct.lenscale,optimState);
        else
            % Distance between vector and set of poll vectors
            dist = zeros(size(U,1),size(ui,1));
            for i = 1:size(ui,1)
                dist(:,i) = udist(U,ui(i,:),gpstruct.lenscale,optimState);
            end
            dist = min(dist,[],2);
        end
        [distord,ord] = sort(dist,'ascend');
        
        % Keep only points within a certain (rescaled) radius from target
        nmin = options.Ndata - options.BufferNdata;
        nmax = min(options.Ndata,optimState.Xmax); 
        radius = 4*gpstruct.effectiveradius;
        nclose = sum(distord <= radius^2);
        
        % [nmin nmax nclose optimState.Xmax]
        
        if optimState.Xmax <= nmin || ((nclose - nmin) <= options.BufferNdata)
            index = 1:nmax;
            gpstruct.x = U(ord(index),:);
            gpstruct.y = Y(ord(index),:);            
        else            
            index = 1:nmin;
            gpstruct.x = U(ord(index),:);
            gpstruct.y = Y(ord(index),:);        
                                    
            index = size(gpstruct.x,1)+1:nclose;            
            xc = U(ord(index),:);
            yc = Y(ord(index),:);
                        
            if 0
                %idx = kmeans([xc,yc],options.BufferNdata);
                % Periodic variables are not fully supported for clustering
                idx = kmeans(bsxfun(@rdivide,xc,gpstruct.lenscale),options.BufferNdata);
                nClusters = max(idx);
                for iCluster = 1:nClusters
                    subset = (idx == iCluster);
                    xx = xc(subset,:);
                    yy = yc(subset);
                    if 0
                        [~,pos] = min(yy);
                    else
                        pos = randi(numel(yy));
                    end
                    gpstruct.y(end+1) = yy(pos);
                    gpstruct.x(end+1,:) = xx(pos,:); 
                end
            else
                idx = randperm(size(xc,1));
                gpstruct.x = [gpstruct.x; xc(idx(1:options.BufferNdata),:)];
                gpstruct.y = [gpstruct.y; yc(idx(1:options.BufferNdata))];
            end
            
        end
        
    case 'grid'
                
        % GP-based vector scaling
        vv = eye(D);
        %vv = grammSchmidt(randn(D,1))';
        vv = [vv; -vv];        
        % vv = [vv; 0.5*vv; 0.25*vv; -0.25*vv; -0.5*vv; -vv];
        
        if options.HessianUpdate
            % vv = vv;
            Bnew = vv*optimState.C';
            % Global vector normalization
            M = sqrt(sum(Bnew(1:D,:).*Bnew(1:D,:),2));
            N = exp(log(M) - mean(log(M)))./M;
            Bnew = bsxfun(@times, Bnew, repmat(N,[size(Bnew,1)/D,1]));
            vv = bsxfun(@times,Bnew*optimState.meshsize,optimState.scale);
            % vv = [vv; bsxfun(@times,[eye(D);-eye(D)]*optimState.meshsize,optimState.scale.*gpstruct.pollscale)];
        else
            vv = bsxfun(@times,vv*optimState.meshsize,optimState.scale.*gpstruct.pollscale);
        end
        %xc = origunits(uc,optimState);
        %xg = periodCheck(bsxfun(@plus,xc,vv),optimState.LB,optimState.UB,optimState,0);
        %xg = xCheck(xg,optimState.LB,optimState.UB,options.TolMesh,optimState,0);
        %ug = [uc; gridunits(xg,optimState)];

        ug = periodCheck(bsxfun(@plus,uc,vv),optimState.LB,optimState.UB,optimState);
        ug = uCheck(ug,options.TolMesh,optimState,1);
        ug = [uc; ug];
                
        % Distance between vector and set of poll vectors
        dist = zeros(size(U,1),size(ug,1));
        for i = 1:size(ug,1)
            dist(:,i) = udist(U,ug(i,:),gpstruct.lenscale,optimState);
        end
        [~,closest] = min(dist,[],1);        
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');

        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*gpstruct.effectiveradius;
        ntrain = min(options.Ndata, sum(distord <= radius^2));
        
        % Minimum number of points to keep
        ntrain = max([options.MinNdata,options.Ndata-options.BufferNdata,ntrain]);
        
        % Up to the maximum number of available points
        ntrain = min(ntrain, optimState.Xmax);
                
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Cluster observations
        index = 1:ntrain;
        
        % Add closest point
        extraidx = find(any(bsxfun(@eq, ord, closest),2));
        index = [index, extraidx'];
        
        % Add safeguarded points
        for d = 1:D
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
            
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end
        
    case 'covgrid'
        
        if ~options.HessianUpdate
            error('HESSIANUPDATE must be on.');
        end
        
        % GP-based vector scaling
        vv = eye(D);
        vv = [vv; -vv];
        
        Bnew = vv*optimState.C';
        % Global vector normalization
        M = sqrt(sum(Bnew(1:D,:).*Bnew(1:D,:),2));
        N = exp(log(M) - mean(log(M)))./M;
        Bnew = bsxfun(@times, Bnew, repmat(N,[size(Bnew,1)/D,1]));
        vv = bsxfun(@times,Bnew*optimState.meshsize,optimState.scale);

        %xc = origunits(uc,optimState);
        %xg = periodCheck(bsxfun(@plus,xc,vv),optimState.LB,optimState.UB,optimState,0);
        %xg = xCheck(xg,optimState.LB,optimState.UB,options.TolMesh,optimState,0);
        %ug = [uc; gridunits(xg,optimState)];

        ug = periodCheck(bsxfun(@plus,uc,vv),optimState.LB,optimState.UB,optimState);
        ug = uCheck(ug,options.TolMesh,optimState,1);
        ug = [uc; ug];        
        
        % Distance between vector and set of poll vectors
        dist = zeros(size(U,1),size(ug,1));
        for i = 1:size(ug,1)
            dist(:,i) = ugdist(U,ug(i,:),optimState.B,optimState);
        end
        [~,closest] = min(dist,[],1);        
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');
        
        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius;
        ntrain = max(options.MinNdata,max(options.Ndata-options.BufferNdata, min(options.Ndata, sum(distord <= radius^2))));
        ntrain = min(ntrain, optimState.Xmax);
        % ntrain = max(floor(ntrain/2),min(ntrain, sum(distord <= radius^2)));
        
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Cluster observations
        index = 1:ntrain;
        
        % Add closest point
        extraidx = find(any(bsxfun(@eq, ord, closest),2));
        index = [index, extraidx'];
        
        % Add safeguarded points
        for d = 1:D
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
            
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);           
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end

    case 'nogrid'
                
        ug = uc;
                
        % Distance between vector and set of poll vectors
        dist = zeros(size(U,1),size(ug,1));
        for i = 1:size(ug,1)
            dist(:,i) = udist(U,ug(i,:),gpstruct.lenscale,optimState);
        end
        [~,closest] = min(dist,[],1);        
        dist = min(dist,[],2);
        [distord,ord] = sort(dist,'ascend');

        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*gpstruct.effectiveradius;
        ntrain = min(options.Ndata, sum(distord <= radius^2));
        
        % Minimum number of points to keep
        ntrain = max([options.MinNdata,options.Ndata-options.BufferNdata,ntrain]);
        
        % Up to the maximum number of available points
        ntrain = min(ntrain, optimState.Xmax);
                
        % sqrt(distord([ntrain-1,ntrain,min(numel(distord),ntrain+1)])')
        
        % Cluster observations
        index = 1:ntrain;
                
        % Add safeguarded points
        for d = 1:D
            idx1 = find(U(ord, d) < uc(d), 1);
            idx2 = find(U(ord, d) > uc(d), 1);
            index = [index, idx1, idx2];
        end
            
        index = unique(index);
        nextra = numel(index) - ntrain;
        % if nextra > 0; nextra, end        
        
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end
        
        
    case 'neighborhood'
        
        if ~options.HessianUpdate
            error('HESSIANUPDATE must be on.');
        end
                        
        dist = ugdist(U,uc,optimState.B,optimState);
        [distord,ord] = sort(dist,'ascend');
        
        % Keep only points within a certain (rescaled) radius from target
        radius = options.gpRadius*optimState.meshsize;
        ntrain = max(options.MinNdata,min(options.Ndata, sum(distord <= radius^2)));
        ntrain = min(ntrain, optimState.Xmax);
        
        %distord
        %radius
        
        index = 1:ntrain;                    
        gpstruct.x = U(ord(index),:);
        gpstruct.y = Y(ord(index),:);           
        if isfield(optimState,'S'); gpstruct.sd = S(ord(index),:); end        
        
    case 'global'

        globalNdata = min(options.globalNdata, optimState.Xmax);
        
        if optimState.Xmax > globalNdata                
            gpstruct.x = [];
            gpstruct.y = [];
            % Periodic variables are not fully supported for clustering
            idx = kmeans(bsxfun(@rdivide,U,gpstruct.lenscale),globalNdata);
            nClusters = max(idx);
            for iCluster = 1:nClusters
                subset = (idx == iCluster);
                xx = U(subset,:);
                yy = Y(subset);
                if 0
                    gpstruct.y(end+1) = mean(yy);
                    gpstruct.x(end+1,:) = mean(xx,1); 
                else
                    if 1
                        [~,pos] = min(yy);
                    else
                        pos = randi(numel(yy));
                    end
                    gpstruct.y(end+1) = yy(pos);
                    gpstruct.x(end+1,:) = xx(pos,:); 
                end
            end
        else
            gpstruct.x = U;
            gpstruct.y = Y;            
        end
        
        gpstruct.y = gpstruct.y(:);
end

% Transformation of objective function
if options.FitnessShaping
    [gpstruct.y,gpstruct.nonlinf,gpstruct.nonlinmu,gpstruct.deltay] = ...
        fitnessTransform(gpstruct.y);
end

% Check for infinities
erry = ~isfinite(gpstruct.y);
if any(erry)
    ygood = gpstruct.y(~erry);
    ypenalty = max(ygood);
    gpstruct.y(erry) = ypenalty;
end

if strcmpi(method,'add')
    gpstruct.erry = gpstruct.erry | erry;
else
    gpstruct.erry = erry;
end


% Store test points
gpstruct.xi = ui;

% Update gp definitions
gplik = [];
gpstruct.x0 = [];
if rotate_gp && ~strcmpi(method,'add') && ~strcmpi(method,'neighborhood')  % Rotate dataset
    if retrain
        if isfield(gpstruct,'Cinv')
            Cinvold = gpstruct.Cinv;
        else
            Cinvold = eye(D);
        end
        C = inv(optimState.C)*optimState.meshsize;
        gpstruct.C = C;
        gpstruct.Cinv = optimState.C/optimState.meshsize;
        idx = gpstruct.ncovoffset+(1:D);
        covlen = exp(gpstruct.hyp.cov(idx));
        covlen = log(gpstruct.C*Cinvold*covlen) % This needs to be fixed
        for i = 1:D
            covlen(i) = min(max(covlen(i), gpstruct.bounds.cov{i}(1)),gpstruct.bounds.cov{i}(2));
        end
        gpstruct.hyp.cov(idx) = covlen;
        
        % C
    end
    gpstruct.x = gpstruct.x*gpstruct.C';
end
gpstruct = feval(options.gpdefFcn{:},D,gplik,optimState,options,gpstruct,options.gpMarginalize);

% Re-fit Gaussian process (optimize or sample)
if retrain
    
    gpstruct = gpfit(gpstruct,options.gpSamples,options);
    
    % Gaussian process length scale
    if gpstruct.ncovlen > 1
        gpstruct.lenscale = zeros(1,D);
        for i = 1:numel(gpstruct.hyp)
            gpstruct.lenscale = gpstruct.lenscale + gpstruct.hypweight(i)*exp(gpstruct.hyp.cov(gpstruct.ncovoffset+(1:D)))';
        end
        gpstruct.lenscale = gpstruct.lenscale/numel(gpstruct.hyp);
    else
        gpstruct.lenscale = 1;
    end
    
    % gp-based geometric length scale
    ll = options.gpRescalePoll*gpstruct.hyp.cov((1:D)+gpstruct.ncovoffset);
    ll = exp(ll - mean(ll))';
    ll = min(max(ll, optimState.searchmeshsize), (optimState.UB-optimState.LB)./optimState.scale); 
    gpstruct.pollscale = ll;
    
    % gp effective covariance length scale radius
    if options.UseEffectiveRadius
        switch lower(gpstruct.covtype)
            case 'rq'
                alpha = exp(gpstruct.hyp.cov(end));
                gpstruct.effectiveradius = sqrt(alpha*(exp(1/alpha)-1));                
            case 'matern1'
                gpstruct.effectiveradius = 1/sqrt(2);                                   
            case 'matern3'
                % gpstruct.effectiveradius = fzero(@(x)(1+sqrt(3)*x)*exp(-sqrt(3)*x)-exp(-1),1)/sqrt(2);
                gpstruct.effectiveradius = 0.876179713323485;               
            case 'matern5'
                % gpstruct.effectiveradius = fzero(@(x)(1+sqrt(5)*x+5/3*x^2)*exp(-sqrt(5)*x)-exp(-1),1)/sqrt(2);
                gpstruct.effectiveradius = 0.918524648109253;               
            otherwise
                gpstruct.effectiveradius = 1;
        end
    end
    
    % Gaussian process signal variability
    gpstruct.sf = 0;
    for i = 1:numel(gpstruct.hyp)
        gpstruct.sf = gpstruct.sf + gpstruct.hypweight(i)*exp(gpstruct.hyp.cov(gpstruct.ncovoffset+D+1));
    end
    gpstruct.sf = gpstruct.sf/numel(gpstruct.hyp);
    
    % [std(gpstruct.y) gpstruct.sf gpstruct.hyp.lik(1) exp(gpstruct.hyp.lik(2:end)')]
end

try
    % Recompute posterior
    [~,~,~,~,~,gpstruct.post] = mygp(gpstruct.hyp,gpstruct.inf,gpstruct.mean,gpstruct.cov,gpstruct.lik,gpstruct.x,gpstruct.y,uc(1,:));
catch
    % Posterior update failed
    gpstruct.post = [];
end

% gpstruct.hyp.cov(:)'
    
end

%--------------------------------------------------------------------------
function gpstruct = gpfit(gpstruct,Nsamples,options)
%GPFIT Fit Gaussian Process hyper-parameters (optimize or sample).

if isfield(gpstruct,'bounds') && ~isempty(gpstruct.bounds)
    bounds = unwrap(gpstruct.bounds);
    lb = bounds(1:2:end-1);
    ub = bounds(2:2:end);
else
    lb = -Inf;
    ub = Inf;
end

% Initial point #1 is old hyperparameter value
hyp0(1) = gpstruct.hyp;

if Nsamples == 0
    % Initial point #2 is avg of random draw from prior and #1
    
    % Check for possible high-noise mode
    if numel(options.NoiseSize) == 1 || ~isfinite(options.NoiseSize(2)); noise = 1; 
    else noise = options.NoiseSize(2); end    
    highNoise = hyp0(1).lik(1) > (log(options.NoiseSize(1)) + 2*noise);
    
    % Check for mean stuck below minimum
    lowMean = hyp0(1).mean(1) < min(gpstruct.y);
    
    % Conditions for performing a second fit
    secondfit = options.DoubleRefit || highNoise || lowMean;
        
    if secondfit
        hrnd = gppriorrnd(gpstruct.prior,gpstruct.hyp);
        hrnd = 0.5*(unwrap(hrnd) + unwrap(gpstruct.hyp));
        if highNoise; hrnd(end-1) = randn()-2; end   % Retry with low noise magnitude
        if lowMean; hrnd(end) = median(gpstruct.y); end % Retry with mean from median
        hyp0(2) = rewrap(gpstruct.hyp,min(max(hrnd,lb),ub));
    end
    optoptions = optimset('TolFun',0.1,'TolX',1e-4,'MaxFunEval',150);
    % optoptions.Hessian = 'user-supplied';
    
    % Remove errors
    if isfield(gpstruct,'erry') && sum(gpstruct.erry) > 0 && 0
        gpopt = gpstruct;
        gpopt.x(gpstruct.erry,:) = [];
        gpopt.y(gpstruct.erry) = [];
        size(gpopt.x,1)
    else
        gpopt = gpstruct;
    end        
    
    hyp = gpHyperOptimize(hyp0,gpopt,options.OptimToolbox,optoptions,options.NoiseNudge,options.RemovePointsAfterTries);     
    hypw = 1;
else
    [hyp,hypw] = gpHyperSample(hyp0,gpstruct);
end

gpstruct.hyp = hyp;
gpstruct.hypweights = hypw;         % Update samples weigths

% Laplace approximation at MAP solution
if Nsamples == 0 && gpstruct.marginalize
    try
        [~,~,~,~,~,gpstruct.hypHessian] = feval(gpstruct.inf{:},gpstruct.hyp,gpstruct.mean,gpstruct.cov,gpstruct.lik,gpstruct.x,gpstruct.y);
    catch
        gpstruct.hypHessian = [];
    end
end

%[exp(gpstruct.hyp.cov(:))', exp(gpstruct.hyp.lik)]
%gpstruct.hyp.mean

end