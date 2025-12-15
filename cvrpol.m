function raydata = cvrpol(s1raydata,w,z,apcheck)
%CVRPOL returns polarization ray trace data at image surface in CodeV
%
%   function raydata = cvrpol(s1raydata,z,w,apcheck)
%
%   INPUTS: 
%
%   OUTPUT: failsurf = 0 if ray traces to image, or surface number at failure
%           raydata is a vector of image data: X Y Z L M N OP TRN
%           rayfail = 0 if traces to image, 
%                    -N if blocked by aperture or obstruction at surface N
%                    +N if TIR or diffracted into evanescent order at surface N
%
%   See also: CVRGRIDPOL CVR CVRGRID

global CodeV

if nargin<1, s1raydata = [0 0 0 0]'; end % x=y=L/N=M/N= 0.0  (dir tangents)
if nargin<2, w=1; end
if nargin<3, z=1; end
if nargin<4, apcheck=0; end

cvresults = CodeV.RAYPOL(z,w,apcheck,s1raydata);
raydata = cell2mat(cvresults.Item('Output'));

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     