function [LB,UB,PLB,PUB,fixidx] = boundscheck(x0,LB,UB,PLB,PUB)
%BOUNDSCHECK Initial check of bounds and see if there are fixed variables.

nvars = numel(x0);

% Expand scalar inputs to full vectors
if isscalar(LB); LB = LB*ones(1,nvars); end
if isscalar(UB); UB = UB*ones(1,nvars); end
if isscalar(PLB); PLB = PLB*ones(1,nvars); end
if isscalar(PUB); PUB = PUB*ones(1,nvars); end

if isempty(PLB) || isempty(PUB)
    warning('bads:pbUnspecified', 'Plausible lower/upper bounds PLB and/or PUB not specified. Using hard upper/lower bounds instead.');
    if isempty(PLB); PLB = LB; end
    if isempty(PUB); PUB = UB; end
end

% Test that all vectors have the same length
if any([numel(LB),numel(UB),numel(PLB),numel(PUB)] ~= nvars)
    error('All input vectors (X0, LB, UB, PLB, PUB), if specified, need to have the same size.');
end

% Test that all vectors are row vectors
if ~isvector(x0) || ~isvector(LB) || ~isvector(UB) || ~isvector(PLB) || ~isvector(PUB) ...
        || size(x0,1) ~= 1 || size(LB,1) ~= 1 || size(UB,1) ~= 1 || size(PLB,1) ~= 1 || size(PUB,1) ~= 1
    error('All input vectors should be row vectors.');
end

% Test that plausible bounds are finite
if ~all(isfinite(([PLB, PUB]))) 
    error('Plausible interval bounds PLB and PUB need to be finite.');
end

% Test that all vectors are real-valued
if ~isreal([x0; LB; UB; PLB; PUB])
    error('All input vectors should be real-valued.');
end

% Fixed variables (all bounds equal)
fixidx = (x0 == LB) & (LB == UB) & (UB == PLB) & (PLB == PUB);