function surfacelabels = cvsl(s,withnumbers_yn)
%CVSL returns surface number and label in a string
%
%   function surfacelabels = cvsl(s,withnumbers_yn);
%
%   INPUT: s = vector of desired surfaces, default = 1:image
%          withnumbers_yn = logic to add surface number prior to label
%                            default = 1 (yes)
%
%   OUTPUT: surfacelabels = character array of surface label data
%
%   See also CVSD

if nargin<1, s = 1:cvnum; end
if nargin<2, withnumbers_yn = 1; end % default to put numbers in front of label data

surfacelabels = [];
for i=1:length(s)
    label = cveva(['(sll s' int2str(s(i)) ')'],1); %returns a string
        if length(label)<1, label = '-'; end
    if withnumbers_yn
        surfacelabels = strvcat(surfacelabels,[int2str(s(i)) ' ' label]);
    else
        surfacelabels = strvcat(surfacelabels,label);
    end
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     