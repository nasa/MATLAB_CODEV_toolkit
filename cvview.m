function viewdata = cvview(view,s,z,savename)
%CVVIEW draws the current lens in Matlab
%
%   function view = cvview(view,s,z,save);
%
%   INPUTS: view = type of view, 1=YZ, 2=XZ, 3=perspective (default=1)
%           s = surface range to be displayed (default=1:image)
%           z = zoom to plot, (default=1)
%           savename = name of PNG file to save
%
%   OUTPUT: Image of current lens, as a PNG file.
%           If savename is input, the PNG file is saved in the
%           current directory
%
%   See also:

if nargin<1, view = 1; end
if nargin<2, s=1:cvnum; end
if nargin<3, z=1; end
if nargin<4, savename=[]; end

view_type = {'PLC S1 YZ','PLC S1 XZ','VPT S1 -37.8 26.6'};
view = view_type{view};

if ~isempty(savename)
    filename_plt = [savename '.plt'];
else
    filename_plt = 'view.plt';
end
filename_png = [filename_plt(1:end-3) 'png'];

cvcmd(['GRA "' filename_plt '"']); 
cvcmd(['vie; ' view ';' ...
    'sur s' num2str(min(s)) '..' num2str(max(s)) ' z' num2str(z) ';' ...
    ' ;go; GRA T; GCV PNG ' filename_plt]);
cvcmd('GRA T; ');

% pause(1);
viewdata = imread(filename_png);

if nargout<1
    figure('color',[1 1 1]);
    imshow(viewdata, Interpolation="bilinear");
    viewdata = [];
end

delete(filename_plt);
if isempty(savename)
    delete(filename_png);
end
    
% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     