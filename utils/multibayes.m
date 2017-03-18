function [xnew,optimState,gpstruct] = multibayes(x,y,LB,UB,hyp,options)

if nargin < 5; hyp = []; end
if nargin < 6; options = []; end

nvars = size(x,2);

defopts = bps('defaults');

% options.NoiseSize = [1,0.2];    % Use 'large' noise
options.SearchAcqFcn = {@acqLCB, []};
% options.SearchAcqFcn = {@acqNegPI};
options.gpdefFcn = {@gpdefStationaryNew,'matern5',1};

% Initialize algorithm options
options = initoptions(nvars,defopts,options);    

% Initalize mesh
[u0,LB,UB,PLB,PUB,MeshSizeInteger,MeshSize,TolMesh,optimState] = ...
    initmesh(x(end,:),LB,UB,LB,UB,[],options);
options.TolMesh = TolMesh;

% Rescale variables
optimState.U = gridunits(x,optimState);
optimState.Y = y;
optimState.Xmax = size(x,1);

optimState.searchfactor = Inf;
optimState.sdlevel = 3;
optimState.funccount = size(x,1);

if options.FitLik; gplik = []; else gplik = log(options.TolFun); end

% Put all points on a relatively coarse grid
optimState.U = force2grid(optimState.U, optimState);

[U,ia,ic] = unique(optimState.U,'rows');
idx = unique(ic);
Y = zeros(size(U,1),1);
for i = 1:numel(idx)
    Y(idx(i)) = mean(optimState.Y(ic == idx(i)));
end
optimState.U = U;
optimState.Y = Y;
optimState.Xmax = size(optimState.U,1);

% Initialize gp
gpstruct = feval(options.gpdefFcn{:},nvars,gplik,optimState,options,[],options.gpMarginalize);

if ~isempty(hyp); gpstruct.hyp = hyp; end

% Train gp
% Ndata = [50,100,200,500];
Ndata = 500;

for i = 1:numel(Ndata)
    options.globalNdata = Ndata(i);
    gpstruct = gpTrainingSet(gpstruct, ...
        'global', ...
        optimState.U(end,:), ...
        [], ...
        optimState, ...
        options, ...
        1);
end

gpstruct.hyp.cov
gpstruct.hyp.lik
gpstruct.hyp.mean

%[~,~,~,~,~,post] = gp(gpstruct.hyp(1),gpstruct.inf,gpstruct.mean,gpstruct.cov, ...
%    gpstruct.lik,gpstruct.x,gpstruct.y,acqx);
%gpstruct.post = post;

%try
    [~,idx] = min(optimState.Y);
    u0 = optimState.U(idx,:);        
    [unew,optimState] = searchMax(options.SearchAcqFcn,u0,gpstruct,LB,UB,optimState,options);
    
    xnew = origunits(unew,optimState);    
%catch
%    warning('multibayes failed, random initialization.');
%    xnew = rand(1,nvars).*(UB-LB) + LB;
%end

end