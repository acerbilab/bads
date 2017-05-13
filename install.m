% MATLAB installation script for BADS
%
% Copyright (c) by Luigi Acerbi 2017

display('Installing BADS... ');

me = mfilename;                                 % what is my filename
pathstr = fileparts(which(me));                 % get my location
addpath(pathstr);                               % add to the path

try
    savepath;                                   % save path
    display('Installation successful!');
    success_install_flag = 1;
catch
    display('Installation failed: could not save path.');
    display('You need to manually add BADS''s installation folder to your MATLAB search path (and save it).'); 
    display('See the <a href="https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html">MATLAB documentation</a> for more information.'); 
    success_install_flag = 0;
end

if success_install_flag
    type([pathstr filesep 'docs' filesep 'README.txt']);
    fprintf('\n');
end

clear me pathstr success_install_flag