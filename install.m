% MATLAB installation script for BADS
%
% Copyright (c) by Luigi Acerbi 2017

display('Installing BADS... ');

me = mfilename;                                 % what is my filename
pathstr = fileparts(which(me));                 % get my location
addpath(pathstr);                               % add to the path
clear me pathstr
                             
try
    savepath;                                   % save path
    display('Installation successful.\n');
catch
    display('Installation failed: could not save path.');
    display('You need to manually add BADS''s installation folder to your MATLAB search path (and save it).'); 
    display('See the <a href="https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html">MATLAB documentation</a> for more information.'); 
end