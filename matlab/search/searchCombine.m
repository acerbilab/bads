function xs = searchCombine(searchFcns,vec,x,gpstruct,LB,UB,optimState,options)
%SEARCHCOMBINE Combine two or more search methods.

if nargin < 3
    xs = 'com';
    return;
end

% VEC not specified, all arguments shifted by one
if nargin < 8
    options = optimState;
    optimState = UB;
    UB = LB;
    LB = gpstruct;
    gpstruct = x;
    x = vec;
    vec = [];
end    

nF = numel(searchFcns);

% If not specified, divide equally among search functions
if isempty(vec) || ~all(isfinite(vec)); vec = ones(1,nF); end

% Normalize vector
vec = vec/sum(vec);

% Divide number of samples among search functions
Nsearch = diff(round([0,cumsum(vec)]*options.Nsearch));

xs = [];

for i = 1:nF
    options.Nsearch = Nsearch(i);

    if ischar(searchFcns{i}) || isa(searchFcns{i}, 'function_handle');
        searchFcns{i} = {searchFcns{i}};
    end
    
    % Combine candidate search points
    xs = [xs; feval(searchFcns{i}{:},x,gpstruct,LB,UB,optimState,options)];
    
end