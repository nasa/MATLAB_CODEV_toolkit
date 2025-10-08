function map = cvmap(s,z)
%CVMAP  returns max aperture of surface s, zoom z, in CodeV
%
%   function map = cvmap(s,z)
%
%   INPUTS: s = surfaces of interest (default=all surfaces)
%           z = zoom position (default=1)
%
%   OUTPUT: map = vector of maximum apertures for given surfaces
%
%   See also: CVGC CVSD CVSL CVBESTSPH

global CodeV

if nargin<1, s=1:cvnum; end
if nargin<2, z=1; end

for ii=1:length(s)
    map(ii) = invoke(CodeV,'GetMaxAperture',s(ii),z);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     