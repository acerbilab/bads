function v = maskindex(v,idx)
% MASKINDEX Mask non-indexed elements in vector

v(~idx) = 0;
