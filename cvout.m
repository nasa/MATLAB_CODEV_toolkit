function output = cvout()
%CVOUT  retrieves CodeV output
%
%   function  output = cvout()
%
%   Output is text returned from CodeV command line
%
%   See also CVON, CVOFF, CVDB, CVIN, CVOPEN, CVSAVE
%

global CodeV

output = invoke(CodeV,'GetCommandOutput');

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     