function [rows,cols] = cvbufput(data,buffer)
%CVBUFPUT puts array data into the CodeV buffer
%
%   INPUT:  data = array of data
%           buffer = number of buffer, default = b0
%           
%   OUTPUT: rows = 
%           cols = 
%
%   See also CVIN
%  

global CodeV

if nargin<2, buffer=1; end

data = double(data);

[rows,cols] = size(data);
cvcmd(['buf del b' num2str(buffer)]);
[success] = invoke(CodeV, 'ArrayToBuffer', rows, cols, buffer, data);
if success~=0, disp('Data transfer failed! '); return; end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     