function [fval,optimState] = funlogger(fun,u,optimState,state,varargin)
%FUNLOGGER Call objective function and do some bookkeeping.
%   [~,OPTIMSTATE] = FUNLOGGER(FUN,U,OPTIMSTATE,'init') starts logging
%   function FUN with starting point U and optimization struct OPTIMSTATE.
%
%   [~,OPTIMSTATE] = FUNLOGGER(FUN,U,OPTIMSTATE,'init',NMAX) stores up to
%   NMAX function values (default NMAX=1e4).
%
%   [~,OPTIMSTATE] = FUNLOGGER(FUN,U,OPTIMSTATE,'init',NMAX,1) also stores
%   heteroskedastic noise (second output argument) from the logged function.
%
%   [FVAL,OPTIMSTATE] = FUNLOGGER(FUN,U,OPTIMSTATE,'iter') evaluates function
%   FUN at U with optimization struct OPTIMSTATE, return function value FVAL.
%   FUN must take a vector input and returns a scalar value and, optionally,
%   the (estimated) SD of the returned value (if heteroskedastic noise
%   handling is on).
%
%   [FVAL,OPTIMSTATE] = FUNLOGGER(FUN,U,OPTIMSTATE,'single') as 'iter' but
%   does not store function values in the cache.
%
%   [~,OPTIMSTATE] = FUNLOGGER(FUN,U,OPTIMSTATE,'done') finalizes stored
%   function values.

%   Luigi Acerbi 2017

fval = [];
switch lower(state)
    case 'init' % Start new function logging session
    
        nvars = numel(u);

        % Number of stored function values
        if nargin > 4; nmax = varargin{1}; else; nmax = []; end
        if isempty(nmax); nmax = 1e4; end
        
        % Heteroscedastic noise provided as 2nd output argument of FUN?
        if nargin > 5; henoise_flag = varargin{2}; else; henoise_flag = []; end
        if isempty(henoise_flag); henoise_flag = false; end

        optimState.funccount = 0;

        % Passing empty or new OPTIMSTATE
        if ~isfield(optimState,'X') || isempty(optimState.X)    
            optimState.X = NaN(nmax,nvars);
            optimState.U = NaN(nmax,nvars);
            optimState.Y = NaN(nmax,1);
            if henoise_flag; optimState.S = NaN(nmax,1); end
            optimState.Xn = 0;                    % Last filled entry
            optimState.Xmax = 0;                  % Maximum entry index

        else % Receiving previous evaluations (e.g. from previous run)
            
            % Find last evaluated point
            optimState.Xmax = find(~isnan(optimState.Y),1,'last');
            if ~isfield(optimState,'Xn') || isempty(optimState.Xn)
                optimState.Xn = optimState.Xmax;
            end
            
            % Current cache is smaller than previous one, keep most recent
            % evaluated points
            if optimState.Xmax > nmax
                optimState.X = circshift(optimState.X,-optimState.Xn);
                optimState.X = optimState.X(end-nmax+1:end,:);
                optimState.Y = circshift(optimState.Y,-optimState.Xn);
                optimState.Y = optimState.Y(end-nmax+1:nmax,:);
                if henoise_flag
                    optimState.S = circshift(optimState.S,-optimState.Xn);
                    optimState.S = optimState.S(end-nmax+1:nmax,:);
                end
                optimState.Xn = 0;
                optimState.Xmax = nmax;
            else
                offset = nmax - size(optimState.X,1);
                optimState.X = [optimState.X; NaN(offset,nvars)];
                optimState.Y = [optimState.Y; NaN(offset,1)];
                if henoise_flag
                    optimState.S = [optimState.S; NaN(offset,1)];
                end                
            end
            optimState.U = gridunits(optimState.X,optimState);
        end
        optimState.funevaltime = NaN(nmax,1);
        optimState.totalfunevaltime = 0;
    
    case {'iter','single'} % Evaluate function (and store output for 'iter')

        x = origunits(u,optimState);    % Convert back to original space
        % Heteroscedastic noise?
        henoise_flag = isfield(optimState,'S');
        
        try
            t_fun = tic;
            if henoise_flag
                [fval,fsd] = fun(x);
            else
                fval = fun(x);
            end
            t = toc(t_fun);
            
            % Check returned function value
            if ~isscalar(fval) || ~isfinite(fval) || ~isreal(fval)
                error(['The returned function value must be a finite real-valued scalar (returned value: ' mat2str(fval) ').']);
            end
            
            % Check returned function SD
            if henoise_flag && (~isscalar(fsd) || ~isfinite(fsd) || ~isreal(fsd) || fsd <= 0.0)
                error(['The returned estimated SD (second function output) must be a finite, non-negative real-valued scalar (returned SD: ' mat2str(fsd) ').']);
            end            
            
        catch fun_error
            warning('bads:funError', ['Error in executing the logged function ''' func2str(fun) ''' with input: ' mat2str(x)]);
            rethrow(fun_error);
        end
        
        % Update function records
        optimState.funccount = optimState.funccount + 1;
        if strcmpi(state, 'iter')
            optimState.Xn = max(1,mod(optimState.Xn+1, size(optimState.X,1)));
            optimState.Xmax = min(optimState.Xmax+1, size(optimState.X,1));
            optimState.X(optimState.Xn,:) = x;
            optimState.U(optimState.Xn,:) = u;
            optimState.Y(optimState.Xn) = fval;
            if henoise_flag
                optimState.S(optimState.Xn) = fsd;            
            end
            optimState.funevaltime(optimState.Xn) = t;
        end
        optimState.totalfunevaltime = optimState.totalfunevaltime + t;
        
    case 'done' % Finalize stored table
        
        henoise_flag = isfield(optimState,'S');
        
        if optimState.Xmax < size(optimState.X,1)
            optimState.X = optimState.X(1:optimState.Xmax,:);
            optimState.Y = optimState.Y(1:optimState.Xmax);
            if henoise_flag; optimState.S = optimState.S(1:optimState.Xmax); end
            optimState.funevaltime = optimState.funevaltime(1:optimState.Xmax);
        else
            optimState.X = circshift(optimState.X,-optimState.Xn);
            optimState.Y = circshift(optimState.Y,-optimState.Xn);
            if henoise_flag; optimState.S = circshift(optimState.S,-optimState.Xn); end
            optimState.funevaltime = circshift(optimState.funevaltime,-optimState.Xn);        
        end
        optimState = rmfield(optimState,'U');
        
    otherwise        
        error('Unknown FUNLOGGER action.');
end

end