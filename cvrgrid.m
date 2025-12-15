function [data,mask] = cvrgrid(nrd,f,s,w,z)
%CVRGRID uses CODE V RAYTRA engine to calculated a grid of rays.
%
% function [data,mask] = cvrgrid(nrd,f,s,w,z);
%
%   INPUTS: nrd = gridsize, default=32
%           f = field, default = 1
%           s = surface, default = image surface
%           w = wavelength number, default = 1
%           z = zoom number,default = 1
%
%   OUTPUT: data = 3D grid of ray data:  x=1, y=2, z=3, L=4, M=5, N=6,
%                   OPL=7, transmission=8, raytra return flag =9
%           mask = Defines the pupil mask based on a ray failure map.
%           
%   See also, CVRAYGRID 

if nargin<1, nrd = 32; end
if nargin<2, f = 1; end
if nargin<3, s = cvnum; end
if nargin<4, w = 1; end
if nargin<5, z = 1; end

cvmacro('cvrgrid.seq', [nrd,f,s,w,z]);
data = cvbuf(1,1:nrd*9,1:nrd);

data = reshape(data',nrd,nrd,9);
mask = data(:,:,1)==0; %CodeV returns 0 for ray pass
data = data(:,:,[2:end 1]); %put raytra return data last

if nargout<1
    figure,
    imagesc(data(:,:,7).*mask);
    colorbar, axis tight, axis square
    title('OPL')
    data = [];
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     