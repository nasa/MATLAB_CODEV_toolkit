function [CODEVpath,defaultspath,CODEVserver] = cvsetup(newconfigfile)
% Creates a config file, and returns location of CODE V installation and defaults.seq
% function [CODEVpath,defaultspath,CODEVserver] = cvsetup(newsetup)
%   INPUT:  newconfigfile = true, to setup a new config file, default = false (no)
%                         = false, to load data from cvconfigfile.m
%                         = -1, to use hard-coded locations in this function

if nargin<1, newconfigfile = false; end

if newconfigfile==-1 
    %USER defined software version and installation locations should be put here
    CODEVpath = 'C:\CODEV115';
    defaultspath = 'C:\CVUSER\defaults.seq';
    CODEVserver = 'CODEV.Application.115';   
elseif ~newconfigfile % use setup locations in configfile
    fname_config = [prefdir filesep 'cv.mat'];
    if newconfigfile || ~exist(fname_config, 'file') %setup a new configfile if needed
        cvconfigfile;
    end
    load(fname_config,'CODEVpath','defaultspath','CODEVserver');
else
    cvconfigfile; % setup a new config file
    fname_config = [prefdir filesep 'cv.mat'];
    load(fname_config,'CODEVpath','defaultspath','CODEVserver');
end

% Get current version from CODE V only if no output requested
if nargout<1
    disp(['You are currently running CODE V version:  ' CODEVserver ]);
    disp('You are currently running CODE V toolkit:  2021a');
    [~,matlab_date] = version;
    if datenum(matlab_date) < datenum('January 1, 2021')
        disp(['WARNING: CODE V toolkit requires MATLAB release >= 2021a']);
    end   
end

