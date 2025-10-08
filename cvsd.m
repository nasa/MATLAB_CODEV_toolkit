function surfacedata = cvsd(s)
%CVSD returns surface data in a cell array
%   function surfacedata = cvsd(s);
%   INPUT: s = vector of desired surfaces, default = 1:image
%   OUTPUT: surfacedata = cell array of surface data
%
%   See also: CVSL

if nargin<1, s = 1:cvnum; end

surfacedata = {};
for i=1:length(s)
    sstr = int2str( s(i) );
    surfacedata{i,1} = cvsl(s(i));
    surfacedata{i,2} = cveva(['(typ surf s' sstr ')'],1);
    surfacedata{i,3} = cveva(['(rdy s' int2str(s(i)) ')']);
    surfacedata{i,4} = cveva(['(thi s' int2str(s(i)) ')']);
    surfacedata{i,5} = cveva(['(map s' int2str(s(i)) ')']);
    surfacedata{i,6} = cveva(['(gla s' sstr ')'],1);
    surfacedata{i,7} = cveva(['(rmd s' sstr ')'],1);
end

if nargout<1
    sdlabels = {'Label','Type','RDY','THI','MAP','GLA','RMD'};
    surfacedata = cat(1,sdlabels,surfacedata);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     