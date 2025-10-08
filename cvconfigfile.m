function cvconfigfile

% cvconfigfile
% Easier / automatic way to configure the NASA MATLAB to CODE V Toolkit
% than having to remember to go back an edit a couple of files after you
% upgrade CODE V
%
% DISCUSSION:
% Config info is stored in cv.mat in the prefdir folder, and includes:
%       includes:   CODEVpath    = 'C:\CODEV113'
%                   CODEVserver  = 'CODEV.Command'
%                   defaultspath = 'Z:\0docs\CVUSER\defaults.seq'
% 
% HISTORY:
% 1/14/2019 David Aronstein, Initial release
% 2023.06.12 Joe Howard, minor updates

fname_config = [prefdir filesep 'cv.mat'];

if exist(fname_config, 'file')
    load(fname_config);
end

% Find out if user wants to select partcular version of CODEV to use

which_codev = questdlg('How should CODE V be configured?', 'CODE V Config', 'Use Most Recent', 'Select Version to Use', 'Use Most Recent');

%CODEVpath and CODEVserver
if strcmp(which_codev, 'Use Most Recent')

    % Find CODE V verson
    dir_list = dir('C:\CODEV*');
    indx = 1;
    if (numel(dir_list) > 1)
        names = {dir_list(:).name};
        [~, sort_indx] = sort(names);
        indx = sort_indx(end);
    end
    recentversion = dir_list(indx).name;
    CODEVpath = ['C:\' recentversion];
%     CODEVserver = ['CODEV.Application.' recentversion(6:end)] ;
    CODEVserver = 'CODEV.Application'; % starting with CodeV202203
    
else %choose folder
    
    CODEVpath = 0;
    while isnumeric(CODEVpath)
        CODEVpath = uigetdir('C:\', 'Select folder with desired CODE V version');
    end
    % Get the folder name in isolation
    % indx = strfind(CODEVpath, filesep);
    % dir_name = CODEVpath((indx(end)+1):end);
    % CODEVserver = ['CODEV.Command.' strrep(dir_name, 'CODEV', '')];
    % CODEVserver = ['CODEV.Application.' strrep(dir_name, 'CODEV', '')];
    CODEVserver = 'CODEV.Application'; % starting with CodeV202203
    
end

% CVUSER folder and defaultspath
defaultspath_dir = getenv('USERPROFILE');
if exist('defaultspath', 'var')
    indx = strfind(defaultspath, filesep);
    defaultspath = defaultspath(1:(indx(end)-1));
    if exist(defaultspath, 'dir')
        defaultspath_dir = defaultspath;
    end
end
defaultspath = 0;
while isnumeric(defaultspath)
    defaultspath = uigetdir(defaultspath_dir, 'Select CODE V CVUSER folder');
end
defaultspath = [defaultspath filesep 'defaults.seq'];
if ~exist(defaultspath,"file")
    disp('WARNING: defaults.seq file not found in CVUSER directory...');
    [filename,pathname] = uigetfile({'*.seq', 'SEQ Files (*.seq)';...
                                     '*.*', 'All files (*.*)'},'Choose a defaults file'); 
    defaultspath = [pathname filename];
end

% Save cv.mat config file
save(fname_config, 'CODEVserver', 'CODEVpath', 'defaultspath');
