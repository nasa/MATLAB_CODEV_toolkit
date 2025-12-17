function cvreadme()
% CVREADME gives a basic introduction to the Matlab CODE V toolkit
% 
% INSTALLATION NOTES:
% 1. CODE V must be installed on the computer to use this toolkit.
% 2. cvsetup.m must be modified to reflect the version of CodeV installed
%   (e.g. version 10.2)
% 3. While multiple sessions of Code V can be called, care should be used
%    in running parallel sessions since only a single global CodeV variable is
%    expected.
% 
% USAGE NOTES:
% 1. Code V runs in the background, so you will not see a Code V GUI window in Windows.
% 2. If Code V or Matlab doesn't respond, open the task manager and kill the "CVCOMM~2.EXE" 
%   process.  This will shut down Code V, but should not affect your Matlab session.
% 3. Use function "cvhelp" to list Code V toolkit functions.
%
% EXAMPLE SESSION: (type the following commands in order at the matlab
% command prompt)
% >>cvon                % starts Code V, puts global "CodeV" variable in workspace
% >>cvhelp              % lists functions and help for CodeV toolkit
% >>cvin                % opens a gui to load a .seq file into Code V
% >>opd0 = cvpma;       % gets pupil map data from Code V
% >>cvshift(2,1,1);     % moves surface 2 in x direction by 1 lens unit
% >>opd1 = cvpma;       % gets new pupil map data from Code V
% >>cvoff               % kills CodeV session
% 
% 
% Please email comments and questions to Joseph.M.Howard@nasa.gov
% 
% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2008.03.19     

help cvreadme

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2010.10.05    