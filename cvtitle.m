function title = cvtitle(z)
%
%CVTITLE queries titles in CODE V
%
%   function title = cvtitle(z)
%
%   INPUT:  z = desired titles zoom positions, default = all zooms
%   OUTPUT: titles = names of zoom positions
%
%   See also CVSD, CVF, CVW, CVZ
%

if nargin<1, z = 1; end

[nums,numf,numw,numz] = cvnum;

for i=1:numz
    title{i,1} = cveva(['(title z' int2str(i) ')'],1); 
end

if nargout<1
    title_table = table((1:numz)',title,'VariableNames',{'Zoom Position','Title'});
   disp(title_table);
end

title = title{z};

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     