function [output,cmd] = cvcmd(cmd)
%CVCMD  sends command to CODE V command line
%
%   function  [output,cmd] = cvcmd(cmd)
%
%   cmd = CodeV command line script to be executed
%   Output is text returned from CodeV command line
%
%   See also CVON, CVOFF, CVDB, CVIN, CVOPEN, CVSAVE
%

global CodeV

%          invoke(CodeV,'Command',''); % Send null string to reset CV output
output = invoke(CodeV,'Command',cmd);

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     