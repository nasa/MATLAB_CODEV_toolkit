function cvstop()
%CVSTOP stops a running CodeV command
%
%   function cvstop()
%
%   See also CVON, CVOFF, CVDB, CVIN, CVOPEN, CVSAVE
%

global CodeV
CodeV.StopCommand;

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     