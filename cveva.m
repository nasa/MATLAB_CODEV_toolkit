function output = cveva(cmd,output_type)
%CVEVA sends command to CODE V command line over existing COM link
%       Output is number converted from string return from CodeV command line
%       NOTE: max precision of output is limited due to string conversion,
%             if float precision is required, use CVEVAF function
%
%   function output = cveva(cmd,output_type)
%
%   cmd = CODE V command line script to be executed in single quotes: 'x'
%   output_type: 0=number (default), 1=string
% 
%   See also CVCMD, CVDB, CVSD, CVEVA
%

global CodeV

if nargin < 2,
    output_type = 0;
end

if output_type == 0, 
    cvcmd(['buf put b1 i1 j1 ' cmd '; buf put b1 i1 j2 0.0;']); %put data into b1 i1 j1
    data = cvbuf(1,1,1:2); %cvbuf only works in getting more than one buffer element 
    output = data(1);

elseif output_type == 1, 
    output = invoke(CodeV,'EvaluateExpression',cmd);
    
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     