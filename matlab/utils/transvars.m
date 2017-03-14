function varargout = transvars(varargin)
%TRANSVARS Standardize variables via linear or nonlinear transformation
%
%  TRINFO = TRANSVARS(NVARS,LB,UB,PLB,PUB) returns the transformation 
%  structure TRINFO for a problem with NVARS dimensions, lower/upper bounds
%  respectively LB and UB and plausible lower/upper bounds respectively PLB 
%  and PUB. LB and UB are either scalars or row arrays that can contain 
%  both real numbers and Inf's. PLB and PUB are either scalars or row arrays
%  and need to be finite.
%  The ordering LB <= PLB < PUB <= UB needs to hold coordinate-wise.
%
%  The standardized transform maps PLB and PUB to the hypercube [-1,1]^NVARS.
%  If PLB and/or PUB are empty, LB and/or UB are used instead. Note that
%  at least one among LB, PLB and one among UB, PUB needs to be nonempty.
%
%  Y = TRANSVARS(X,'dir',TRINFO) performs direct transform of array X into
%  standardized array Y according to transformation encoded in structure 
%  TRINFO. X must be a N x NVARS array, where N is the number of input data 
%  and NVARS is the number of dimensions.
%
%  X = TRANSVARS(Y,'inv',TRINFO) performs the inverse transformation from
%  standardized array Y to original array X.
%
%
%  Author: Luigi Acerbi
%  e-mail: luigi.acerbi@gmail.com
%  Release: 0.9
%  Release date: 15/Jan/2016

NumEps = 1e-6;         % Accepted numerical error
NumMax = sqrt(realmax);     % Maximum encoded number

% MaxPrecision = 17; % Maximum precision for a double

if nargin < 3
    error('TRANSVARS requires a minimum of three input arguments.');
end

%% Transform variables
if isstruct(varargin{3})

    trinfo = varargin{3};    

    if isempty(trinfo)
        varargout{1} = varargin{1}; % Return untransformed input
    else
        
        direction = varargin{2};
        if isempty(direction)
            error('The transformation direction cannot be empty. Allowed values are direct (''dir'' or ''d'') and inverse (''inv'' or ''i'').');
        end
        
        if lower(direction(1)) == 'd'
            x = varargin{1};
            if ischar(trinfo.g); g = str2func(trinfo.g); else g = trinfo.g; end
            y = g(x);
            y = min(max(y,trinfo.lb),trinfo.ub);    % Force to stay within bounds

            varargout{1} = reshape(y,size(x));
            
        elseif lower(direction(1)) == 'i'
            y = varargin{1};
            if ischar(trinfo.ginv); ginv = str2func(trinfo.ginv); else ginv = trinfo.ginv; end
            x = ginv(y);
            x = min(max(x,trinfo.oldbounds.lb),trinfo.oldbounds.ub);    % Force to stay within bounds
            varargout{1} = reshape(x,size(y));
            
        else
            error(['Unkwnown transformation direction ''' direction '''. Allowed values are direct (''dir'' or ''d'') and inverse (''inv'' or ''i'').']);
        end
    end    
    
else
%% Create transform

    nvars = varargin{1};
    lb = varargin{2}(:)';
    ub = varargin{3}(:)';
    if nargin > 3; plb = varargin{4}(:)'; else plb = []; end
    if nargin > 4; pub = varargin{5}(:)'; else pub = []; end
    if nargin < 6; trinfo.logct = []; else trinfo.logct = varargin{6}(:)'; end
    
    % Plausible bounds are equal to hard bounds if not specified
    if isempty(plb); plb = lb; end
    if isempty(pub); pub = ub; end
    
    if isempty(lb) && isempty(plb) || isempty(ub) && isempty(pub)
        error('At least one among LB, PLB and one among UB, PUB needs to be nonempty.');
    end
    
    % Empty LB and UB are Infs
    if isempty(lb); lb = -Inf; end
    if isempty(ub); ub = Inf; end
    
    % Convert scalar inputs to row vectors
    if isscalar(lb); lb = lb*ones(1,nvars); end
    if isscalar(ub); ub = ub*ones(1,nvars); end
    if isscalar(plb); plb = plb*ones(1,nvars); end
    if isscalar(pub); pub = pub*ones(1,nvars); end
    
    % Check finiteness of plausible range
    assert(all(isfinite(([plb, pub]))), ...
        'Plausible interval ranges PLB and PUB need to be finite.');

    % Check that the order of bounds is respected
    assert(all(lb <= plb & plb < pub & pub <= ub), ...
        'Interval bounds needs to respect the order LB <= PLB < PUB <= UB for all coordinates.');
    
    % Nonlinear log transform
    if isempty(trinfo.logct)
        trinfo.logct = NaN(1,nvars);
    elseif isscalar(trinfo.logct)
        trinfo.logct = trinfo.logct*ones(1,nvars);
    end
    
    % A variable is converted to log scale if all bounds are positive and 
    % the plausible range spans at least one order of magnitude
    for i = find(isnan(trinfo.logct))
        trinfo.logct(i) = all([lb(i) ub(i) plb(i) pub(i)] > 0) & (pub(i)/plb(i) >= 10);    
    end    
    trinfo.logct = logical(trinfo.logct);
    
    % Convert infinities to very large numbers
    
    % Transform to log coordinates
    trinfo.oldbounds.lb = lb;
    trinfo.oldbounds.ub = ub;
    trinfo.oldbounds.plb = plb;
    trinfo.oldbounds.pub = pub;    
    trinfo.lb = lb; trinfo.ub = ub; trinfo.plb = plb; trinfo.pub = pub;
    trinfo.lb(trinfo.logct) = log(trinfo.lb(trinfo.logct));
    trinfo.ub(trinfo.logct) = log(trinfo.ub(trinfo.logct));
    trinfo.plb(trinfo.logct) = log(trinfo.plb(trinfo.logct));
    trinfo.pub(trinfo.logct) = log(trinfo.pub(trinfo.logct));

    trinfo.mu = 0.5*(trinfo.plb + trinfo.pub);
    trinfo.gamma = 0.5*(trinfo.pub - trinfo.plb);
        
    z = ['maskindex( (x - ' vec2str(trinfo.mu) ') ./ ' vec2str(trinfo.gamma) ' , ' vec2str(~trinfo.logct) ')'];
    zlog = ['maskindex( ( (log(abs(x) + (x == 0)) - ' vec2str(trinfo.mu) ') ./ ' vec2str(trinfo.gamma) ' ), ' vec2str(trinfo.logct) ')'];
    
    switch sum(trinfo.logct)
        case 0
            trinfo.g = ['@(x) (' z ')'];
            trinfo.ginv = ['@(y) ' vec2str(trinfo.gamma) ' .* (y) + ' vec2str(trinfo.mu) ];        
            
        case nvars
            trinfo.g = ['@(x) (' zlog ')'];
            trinfo.ginv = ['@(y) min(realmax,exp(' vec2str(trinfo.gamma) ' .* (y) + ' vec2str(trinfo.mu) '))'];        
            
        otherwise
            %trinfo.g = ['@(x) (1-' vec2str(trinfo.logct) ').*(' z ')' ... 
            %    '+ ' vec2str(trinfo.logct) '.*(' zlog ')'];
            %trinfo.ginv = ['@(y) (1-' vec2str(trinfo.logct) ') .* (' vec2str(trinfo.gamma) ' .* (y) + ' vec2str(trinfo.mu) ') + ' ...
            %    vec2str(trinfo.logct) ' .* min(realmax,exp(' vec2str(trinfo.gamma) ' .* (y) + ' vec2str(trinfo.mu) '))'];        
            trinfo.g = ['@(x) ' z ' + ' zlog];
            trinfo.ginv = ['@(y) maskindex(' vec2str(trinfo.gamma) ' .* (y) + ' vec2str(trinfo.mu) ', ' vec2str(~trinfo.logct) ') + ' ...
                ' maskindex( min(realmax,exp(' vec2str(trinfo.gamma) ' .* (y) + ' vec2str(trinfo.mu) ')), ' vec2str(trinfo.logct) ')'];        
    end
    
    % Convert starting values to transformed coordinates
    trinfo.g = str2func(trinfo.g);
    trinfo.ginv = str2func(trinfo.ginv);
    
    if ischar(trinfo.g); g = str2func(trinfo.g); else g = trinfo.g; end
    if ischar(trinfo.ginv); ginv = str2func(trinfo.ginv); else ginv = trinfo.ginv; end
    
    % Check that the transform works correctly in the range
    lbtest = lb;
    lbtest(~isfinite(lb)) = -1/sqrt(eps);
    ubtest = ub;
    ubtest(~isfinite(ub)) = 1/sqrt(eps);
    ubtest(~isfinite(ub) & trinfo.logct) = 1e6;
    
    t(1) = all(abs(ginv(g(lbtest)) - lbtest) < NumEps);
    t(2) = all(abs(ginv(g(ubtest)) - ubtest) < NumEps);
    t(3) = all(abs(ginv(g(plb)) - plb) < NumEps);
    t(4) = all(abs(ginv(g(pub)) - pub) < NumEps);
    assert(all(t), 'Cannot invert the transform to obtain the identity at the provided boundaries.');
    
    trinfo.lb = g(lb);
    trinfo.ub = g(ub);
    trinfo.plb = g(plb);
    trinfo.pub = g(pub);
    
    varargout{1} = trinfo;
    
end

%--------------------------------------------------------------------------
function s = vec2str(v,idx)
% Convert numerical vector to string

if nargin < 2 || isempty(idx); idx = true(size(v)); end
v(~idx) = 0;    % Ignore unwanted entries

MaxPrecision = 17;  % Maximum precision for a double
if size(v,1) > 1; transp = ''''; else transp = []; end
s = '[';
for i = 1:length(v)-1; s = [s, num2str(v(i),MaxPrecision), ',']; end
s = [s, num2str(v(end),MaxPrecision), ']' transp];

% s = ['[' num2str(v(:)',MaxPrecision) ']' transp];



