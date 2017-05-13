function Uprime = uWarp(U,a,b,LB,UB,action)
%UWARP Input warping of U.

if nargin < 6 || isempty(action); action = 'd'; end

% Normalize
U = bsxfun(@rdivide, bsxfun(@minus, U, LB), UB - LB);

U = min(max(U,0),1);

% Warp inputs
switch lower(action(1))
    case 'd'    % Direct transform
        Uprime = 1 - bsxfun(@power, 1 - bsxfun(@power, U, a), b);
    case 'i'    % Inverse transform
        Uprime = bsxfun(@power, 1 - bsxfun(@power, 1 - U, 1./b), 1./a);
    otherwise
        error('Warping ACTION should be (''d'')irect or (''i'')nverse.');
end

Uprime = min(max(Uprime,0),1);

% Stretch back
Uprime = bsxfun(@plus, bsxfun(@times, Uprime, UB - LB), LB);

end