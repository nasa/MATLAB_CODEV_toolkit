function output = cveva(cmd,output_type)
%CVEVA evaluates the input command in CODE V using "EVA"
%   function output = cveva(cmd,output_type)
%
%   INPUT:  cmd, CODE V command line script to be executed in single quotes: 'x'
%           output_type: 0=number (default), 1=string
% 
%   See also CVCMD, CVDB
%

global CodeV

if nargin < 2, output_type = 0; end

output = CodeV.EvaluateExpression(cmd);
    
if output_type == 0
    output = str2double(output);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     