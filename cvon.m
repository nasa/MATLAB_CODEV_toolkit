function CodeV = cvon(path,filename)
%CVON starts the COM link between Matlab and CODE V
%      and generates the workspace variable 'CodeV'
%   function CodeV = cvon(path)
%       INPUTS:  path = starting directory for CODE V, 
%                       or dir structure pointing to a file to load
%
%   USAGE NOTE:  Version of CODE V and path to defaults.seq should be set
%   in the cvsetup.m function.
% 
%   See also CVOFF, CVSETUP, CVCMD, CVDB, CVIN, CVOPEN, CVSAVE 

if ismac, disp('ERROR: Must be run on Windows!!'); return, end

if nargin<1, path = []; end
if nargin<2, filename = []; end

if isstruct(path)
    filename = path.name;
    path = path.folder;
end
if isempty(path), path = cd; end

if exist('CodeV','var')
    disp('CODE V appears to be running, since the variable ''CodeV'' is in use.');
    disp('This version of the CODE V toolkit runs only one instance of CodeV.'); 
    return
end

global CodeV %put global variable onto workspace

[~,defaultsfile,CODEVserver] = cvsetup; %calls USER DEFINED settings in cvsetup/cvconfigfile
% defaultsfile = 'Z:\0docs\CVUSER\defaults.seq';
% CODEVserver = 'CODEV.Application';

% API version >= 11.3
disp('Starting CODE V...');
CodeV = actxserver(CODEVserver); %start COM link
% CodeV.StartingDirectory = cd; %set current directory to start CODE V
CodeV.StartingDirectory = path; %set current directory to start CODE V
CodeV.StartCodeV; %start session of CODE V
cvcmd(['cd "' path '"']); %bug in version 2022.03

% API version < 11.3
% CODEVpath = 'C:\CODEV113';
% defaultspath = 'Z:\0docs\CVUSER\defaults.seq';
% CODEVserver = 'CODEV.Command';
% invoke(CodeV,'SetStartingDirectory',cd); %set current directory as starting directory
% invoke(CodeV,'StartCodeV'); %start session of CODE V

if exist(defaultsfile, 'file')
    cvin(defaultsfile); %load defaults.seq file
end
cvcmd(['pth seq app "' fileparts(which('cvfw.seq')) '"']); %ensures local installation of CODE V can "see" the macro folder

cvsetup; %verifies that CODE V is running, prints results to command window

if ~isempty(filename)
    cvcmd(['in "' filename '"']);
    cvnum;
end

% clear cmd cvtoolkitversion defaultspath CODEVpath
if nargout<1
    assignin('base',"CodeV",CodeV);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2008.03.19       
