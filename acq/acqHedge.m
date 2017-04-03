function [index,ymu,ys,fmu,fs,fpi,optimState] = acqHedge(xi,target,starget,gpstruct,optimState,options,SufficientImprovement)
%ACQPORTFOLIO Portfolio computation of acquisition functions.

if nargin < 4 || isempty(gpstruct); gpstruct.sf = 1; end
if nargin < 5 || isempty(optimState); optimState.funccount = 0; optimState.sdlevel = 0; end
if nargin < 6 || isempty(options); options = []; end
if nargin < 7 || isempty(SufficientImprovement); SufficientImprovement = 0; end

grad = 0;

% Portfolio of acquisition functions (only EI and LCB)

% sdlevel = [0,1,3,10];
t = optimState.funccount + 1;
sdlevel = optimState.sdlevel;

% Different offsets for EI
% xxi = [0.01,0.1,1]*gpstruct.sf;
% xxi = 0.01*log(optimState.funccount+1);
xxi = 0;

% Different ETA values for LCB
% eta = [0.1,0.2,1];
eta = 0.2;

if nargin < 1
    % Initial vector
    nacq = numel(sdlevel)*numel(xxi) + numel(eta);
    index = zeros(1, nacq);
    index(1) = 10;
    % ymu = {'mpi','mei','lcb'};
    ymu = {'mei','lcb'};
    return;
end

n = size(xi,1);
nvars = size(xi,2);
Nhyp = numel(gpstruct.hyp);

if grad == 1
    error('acqHedge:gradient', ...
        'Gradient of acquisition function is not supported.');
end

%target = target - 0.01*log(t);

if grad
    [ymu,ys2,fmu,fs2,hypw,dymu,dys2,dfmu,dfs2] = gppred(xi,gpstruct,'central');
else
    [ymu,ys2,fmu,fs2,hypw] = gppred(xi,gpstruct);
end
fs = sqrt(fs2);
ys = sqrt(ys2);

index = [];

if any(isnan(fmu)) || any(isnan(fs))
    index = NaN(1,optimState.hedge.n);
    fpi = NaN(size(fmu));
    return;
end

% Probability of improvement and expected improvement
for iSD = 1:numel(sdlevel)
    
    delta = starget*sqrt(sdlevel(iSD)^2 + options.TolFun^2);
    % delta = log(t)*SufficientImprovement;
    
    for iIter = 1:numel(xxi)

        % Probability of improvement
        gammaz = (target - delta - xxi(iIter) - fmu)./fs;
        fpi = 0.5*erfc(-gammaz/sqrt(2));      

        % Expected improvement
        z = -fs.*(gammaz.*fpi + exp(-0.5*(gammaz.^2))/sqrt(2*pi));        

        % Expected squared improvement
        % z = -fs.^2.*((gammaz.^2+1).*fpi + gammaz.*exp(-0.5*(gammaz.^2))/sqrt(2*pi));

        % Average probability of improvement over samples
        try
            fpi = sum(bsxfun(@times,hypw(~isnan(hypw)),fpi(~isnan(hypw),:)),1);
        catch
            fpi = Inf(1,n);
        end    

        % Average EI over samples
        try
            z = sum(bsxfun(@times,hypw(~isnan(hypw)),z(~isnan(hypw),:)),1);
        catch
            z = Inf(1,n);
        end
        
        % [~,index(iIter)] = min(fpi);
        % [~,index(iIter+numel(xxi))] = min(z);
        [~,index(end+1)] = min(-fpi);
        [~,index(end+1)] = min(z);
    end
end

% Lower confidence bound
for iIter = 1:numel(eta)
    delta = 0.1;    
    sqrtbetat = sqrt(eta(iIter)*2*log(nvars*t^2*pi^2/(6*delta)));

    % LCB
    z = fmu - sqrtbetat.*fs;

    % Average LCB over samples
    try
        z = sum(bsxfun(@times,hypw(~isnan(hypw)),z(~isnan(hypw),:)),1);
    catch
        z = Inf(1,n);
    end    
    
    [~,index(end+1)] = min(z);    
end

end