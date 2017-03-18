% MATLAB installation script for BADS
%
% Copyright (c) by Luigi Acerbi 2017

fprintf('Installing BADS... ');

me = mfilename;                                 % what is my filename
pathstr = fileparts(which(me));                 % get my location
addpath(pathstr);                               % add to the path
clear me pathstr
                             
try
    savepath;                                   % save path
    fprintf('Installation successful.\n');
catch
    fprintf('Installation failed: could not save path.\n');
    fprintf('You need to manually add BADS''s ./matlab folder to your MATLAB path.\n');
end