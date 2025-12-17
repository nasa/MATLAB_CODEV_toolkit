function cmd = cvcv(pathfilename)
% CVCV starts CODE V GUI version and loads requested file
% function cmd = cvcv(pathfilename)

if nargin<1, pathfilename = []; end

if ismac, disp('ERROR: Must be run on Windows!!'); return, end;

[CODEVpath,defaultspath] = cvsetup();

if isempty(pathfilename)
    [filename,pathname] = uigetfile({'*.seq', 'SEQ Files (*.seq)';...
                                     '*.len', 'LEN Files (*.len)';...
                                     '*.*', 'All files (*.*)'},'Choose a .seq file to open'); 
    pathfilename = [pathname filename];
elseif isstruct(pathfilename) %dir structure input
    pathfilename = [pathfilename.folder filesep pathfilename.name];
end

if ischar(pathfilename)
    path = fileparts(pathfilename);
    if ~strcmp(path,cd)
        cdnow = cd;
        cd(path);
    end
end

if pathfilename(1) == 0
    cmd = [CODEVpath '\codev &'];
else
    % cmd = [CODEVpath '\codev "' pathfilename '" &'];
    cmd = [CODEVpath '\codev /D="' defaultspath '" "' pathfilename '" &'];
end
dos(cmd);

%add def.seq pointer to my own defaults.seq for convenience
% fid = fopen('def.seq','wt');
% deffilename = ['in "' defaultspath '";'];
% fprintf(fid,'%s\n',deffilename);
% fclose(fid);

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     