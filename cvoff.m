function cvoff
%CVOFF  kills the COM link between CODE V and Matlab
%
%   function cvoff
%
%   See also CVON, CVCMD, CVDB, CVIN, CVOPEN, CVSAVE 
%

global CodeV

invoke(CodeV,'StopCodeV');
evalin('base','clear CodeV');

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22      