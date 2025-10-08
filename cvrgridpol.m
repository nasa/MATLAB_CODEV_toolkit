function [data,mask,success] = cvrgridpol(nrd,f,w,z)
%CVRGRIDPOL returns array of polarization ray data at image surface in CodeV
%
%   function data = cvrgridpol(nrd,f,w,z)
%
%   INPUTS: nrd = gridsize, default=32
%           f = field, default = 1
%           w = wavelength number, default = 1
%           z = zoom number,default = 1
%
%   OUTPUT: data = nrd x nrd x 42 array, data items are X Y Z L M N OPD TRN ...
%                     col 1 = chief ray data
%                     col 2:end = rest of grid
%           rayfail = 1 if chief ray fails, 0 if passes

global CodeV

if nargin<1, nrd = 32; end
if nargin<2, f=1; end
if nargin<3, w=1; end
if nargin<4, z=1; end

numrays = nrd*nrd;
coords = linspace(-1,1,nrd);
[xpupcoords,ypupcoords] = meshgrid(coords);
pupcoords = [xpupcoords(:) ypupcoords(:)]';
raydata = zeros(43,numrays+1);
[success,raydata] = invoke(CodeV,'POLGRID',z,w,f,0,numrays,pupcoords,raydata);

for ii=1:43,
    data(:,:,ii) = reshape(raydata(ii,2:numrays+1),nrd,[]);
end
mask = zeros(nrd);
mask(data(:,:,43)==0)=1;

disp('NOTE: This function cannot currently check apertures!')

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     