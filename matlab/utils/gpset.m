% performs argument checks/transformations similar to those found in
% gp.m from GPML, but Hessian friendly

% Copyright (c) 2014--2015 Roman Garnett.
% Modified 2016 Luigi Acerbi.

% function gpstruct = gpset(gpstruct, hyp, inf, mean, cov, lik, prior, bounds, x, y)
function gpstruct = gpset(gpstruct, varargin)

% Input arguments
args = {'hyp', 'inf', 'mean', 'cov', 'lik', 'prior', 'bounds'};

for idx = 1:numel(varargin)
    % Inputs are x and y
    if isnumeric(varargin{idx})
        if numel(varargin) ~= idx+1 || ~isnumeric(varargin{idx+1})
            error('gpstruct:trainingdata', ...
                'Training data X and Y need to be the last two inputs of GPSET.');
        end        
        gpstruct.x = varargin{idx};
        gpstruct.y = varargin{idx+1};
        break;
    end
    
    if idx > numel(varargin)
        error('gpstruct:inputarguments', ...
            'Too many input arguments.');
    end
    
    % Assign input to struct if different than skip placeholder
    if ~strcmp(varargin{idx}, '~')
        gpstruct.(args{idx}) = varargin{idx};
    end
end

% Default gp values
defaultstruct = struct( ...
    'hyp', struct('cov',[],'lik',[],'mean',[]), ...
    'inf',{@exact_inference}, ...
    'cov', {@isotropic_sqdexp_covariance}, ...
    'lik',{@likGauss}, ...
    'mean', {@zero_mean}, ...
    'x', [], ...
    'y', [], ...
    'prior', [], ...
    'bounds', [] ...
    );

% Assign default values
for field = fields(defaultstruct)'
    if ~isfield(gpstruct,field{:}) || isempty(gpstruct.(field{:}))
        gpstruct.(field{:}) = defaultstruct.(field{:});
    end
end

% Allow string/function handle input; convert to cell arrays if necessary
for field = {'inf', 'cov', 'lik', 'mean'}
    if ischar(gpstruct.(field{:})) || isa(gpstruct.(field{:}), 'function_handle')
        gpstruct.(field{:}) = {gpstruct.(field{:})};
    end
end

% If HYP is a cell array, convert to struct array
if iscell(gpstruct.hyp) && isstruct(gpstruct.hyp{1})
    hyptemp = struct('cov','lik','mean');
    for i = 1:length(gpstruct.hyp)
        hyptemp(i) = gpstruct.hyp{i};
    end
    gpstruct.hyp = hyptemp;
end

% Ensure all hyperparameter fields exist
for field = {'cov', 'lik', 'mean'}
    for i = 1:numel(gpstruct.hyp)
        if (~isfield(gpstruct.hyp(i), field{:}))
            gpstruct.hyp(i).(field{:}) = [];
        end
    end
end

% Check dimension of hyperparameter fields if X is available
D = size(gpstruct.x, 2);
if D > 0
    shortfields = {'cov','lik','mean'};
    longfields = {'covariance','likelihood','mean'};
    for ifield = 1:numel(shortfields)
        field = shortfields{ifield};
        expression = feval(gpstruct.(field){:});
        for i = 1:numel(gpstruct.hyp)
            nparams = numel(gpstruct.hyp(i).(field));
            if (nparams ~= eval(expression))
                error('gpstruct:incorrect_specification', ...
                      'wrong number of %s hyperparameters (%i given, %s expected)', ...
                      longfields{ifield}, ...
                      nparams, ...
                      expression);
            end
        end
    end
end

% Check prior and bounds
for field = {'cov', 'lik', 'mean'}
    % Ensure all prior fields exist
    if (~isfield(gpstruct.prior, field{:}))
        gpstruct.prior.(field{:}) = [];
    end
    % Ensure all bounds fields exist
    if (~isfield(gpstruct.bounds, field{:}))
        gpstruct.bounds.(field{:}) = [];
    end    
    % Assign [-Inf;Inf] bounds by default
    num_entries = numel(gpstruct.hyp(1).(field{:}));
    for i = 1:num_entries
        if numel(gpstruct.prior.(field{:})) < i
            gpstruct.prior.(field{:}){i} = [];
        end
        if numel(gpstruct.bounds.(field{:})) < i || ...
                isempty(gpstruct.bounds.(field{:}){i})
            gpstruct.bounds.(field{:}){i} = [-Inf,Inf];
        end
    end
end

% Add field with (empty) posterior
if ~isfield(gpstruct,'post'); gpstruct.post = []; end
