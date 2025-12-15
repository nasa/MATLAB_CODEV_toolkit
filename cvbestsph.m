function radius = cvbestsph(s,z,ymin,ymax)
%CVBESTSPH  returns radius of best fit sphere
%   function radius = cvbestsph(s,z,ymin,ymax)
%
%   INPUTS: s = surface of interest (default=1)
%           z = zoom position (default=1)
%           ymin = minimum aperture value (default=0)
%           ymax = maximum aperture value (default=map)
%
%

global CodeV

if nargin<1, s=1; end
if nargin<2, z=1; end
if nargin<3, ymin=0.0; end
if nargin<4, ymax=cvmap(s,z); end

curv = CodeV.BESTSPH(s,z,ymin,ymax);
if curv<abs(eps)
    radius = 1e20;
else 
    radius = 1/curv;
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     