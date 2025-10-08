function cmd = cvcv()
% CVCV starts CODE V GUI version and loads file
% function path = cvcv()
% Note that CVUSER directory must be manually placed in this script.

[filename,pathname] = uigetfile({'*.seq', 'SEQ Files (*.seq)';...
                                 '*.len', 'LEN Files (*.len)';...
                                 '*.*', 'All files (*.*)'},'Choose a .seq file to open'); 

cmd = ['C:\CODEV970\codev /DEFAULT=c:\0docs\CVUSER\defaults "' filename '" &'];
dos(cmd);


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     