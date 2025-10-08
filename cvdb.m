function output = cvdb(query,output_type)
%CVDB queries CODE V database item
%     wraps "('cmd')" around the input, and executes cveva
%
%   function  output = cvdb(query,output_type)
%
%   INPUT:  cmd = CODE V command line script to be executed (in quotes)
%           output_type: 0=number (default), 1=string
%   OUTPUT: output = value of database item
% 
%   See also cveva

if nargin<1, help cvdb; return, end
if nargin<2, output_type = 0; end

cmd = ['(' query ')'];
output = cveva(cmd,output_type);

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     