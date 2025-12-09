function cmd = cvcv()
% CVCV starts CODE V GUI version and loads file
% function path = cvcv()

[filename,pathname] = uigetfile({'*.seq', 'SEQ Files (*.seq)';...
                                 '*.len', 'LEN Files (*.len)';...
                                 '*.*', 'All files (*.*)'},'Choose a .seq file to open'); 
pathfilename = [pathname filesep filename];

[CODEVpath,defaultspath] = cvsetup;

cmd = [CODEVpath '\codev /DEFAULT="' defaultspath '" "' pathfilename '" &']; %loads defaults.seq file
%cmd = [CODEVpath '\codev "' filename '" &']; %version with no defaults.seq file

dos(cmd);


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     