% Installation script for BPS
%
% Copyright (c) by Luigi Acerbi Jan-11-2016.

fprintf('Installing BPS... ');

me = mfilename;                                 % what is my filename
pathstr = fileparts(which(me));                 % get my location

dirs = {'acq','gpdef','gpml_fast','init','poll','search','utils'};

addpath(pathstr)
for d = dirs
    addpath([pathstr,filesep,d{:}])
end
clear d dirs me pathstr
                             
try
    savepath;                                   % save path
    fprintf('Installation successful.\n');
catch
    fprintf('Installation failed: could not save path.\n');
    fprintf('You need to manually add all subdirectories to your MATLAB path.\n');
end