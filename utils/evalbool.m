function tf = evalbool(s)
%EVALBOOL Evaluate argument to a bool

if ~ischar(s) % S may not and cannot be empty
        tf = s;
        
else % Evaluation of string S
    if strncmpi(s, 'yes', 3) || strncmpi(s, 'on', 2) ...
        || strncmpi(s, 'true', 4) || strncmp(s, '1 ', 2)
            tf = 1;
    elseif strncmpi(s, 'no', 2) || strncmpi(s, 'off', 3) ...
        || strncmpi(s, 'false', 5) || strncmp(s, '0 ', 2)
            tf = 0;
    else
        try tf = evalin('caller', s); catch
            error(['String value "' s '" cannot be evaluated']);
        end
        try tf ~= 0; catch
            error(['String value "' s '" cannot be evaluated reasonably']);
        end
    end

end

end