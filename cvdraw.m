function cvdraw(view,s,z,savestr)
%CVDRAW draws the current lens in Matlab
%
%   function cvdraw(view,s,z,save);
%
%   INPUTS: view = type of view, 1=YZ, 2=XZ, 3=perspective (default=1)
%           s = surface range to be displayed (default=1:image)
%           z = zoom to plot, (default=1)
%           save = switch to save plot file 0=no, 1=yes, (default=0)
%
%   OUTPUT: Plot of current lens. Opens with the CodeV plot viewer.
%           If the save option is chosen, a plot file is saved in the
%           \graphics subdirectory of the toolkit.
%
%   See also:

if nargin<1, view = 1; end
if nargin<2, s=1:cvnum; end
if nargin<3, z=1; end
if nargin<4, savestr=0; end

view_type = {'PLC S1 YZ','PLC S1 XZ','VPT S1 -37.8 26.6'};
view = view_type{view};

filename = ['view_' datestr(now,'yyyymmdd_HHMMSS') '.plt'];

cvcmd(['GRA ' cvpath '\graphics\' filename ' ; out no']); 
cvcmd(['vie; ' view ';' ...
    'sur s' num2str(min(s)) '..' num2str(max(s)) ' z' num2str(z) ';' ...
    ' ;go;']);
cvcmd('GRA T; out yes');


pause(1);
winopen([ cvpath '\graphics\' filename]);


if savestr==0,
    pause(5);
    delete([ cvpath '\graphics\' filename]);
end
    
% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     