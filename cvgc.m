function    gc = cvgc(s,refsurf,COS)
%CVGC  gets the global coodinates for surfaces, and converts using makeglbs.seq
%
%   function    gc = cvgc(s,refsurf,COS)
%
%   INPUT:  s,z: surfaces and zoom position for coordinates (one zoom only)
%                NOTE: use negative s input to set global coords for s to refsurf 
%           refsurf: enter the reference surface from which the global
%               coordinate system will be specified (def = 1)
%           COS: if true (~=0), orientation is in direction cosines (def = 0)
%                NOTE: not currently implemented in this release.
%           
%   OUTPUT: gc = global coordinates of surface with respect to refsurf
%                X, Y, Z, A(y rot), B(x rot), C(z rot)
%
%   See also: CVSD CVSL CVNUM

ims = cvnum;

if nargin<1, s=1:ims; end
if nargin<2, refsurf =1; end

if sum(s<0)>0, %negative sign sets surfaces s to global coords
    [temp,index] = sort(abs(s(:)),1,'descend');
    for i=1:length(s)
        if s(index(i)) < 0 %set only negative s values
            if -s(index(i)) > refsurf
                cvmacro('makeglbs.seq',[-s(index(i)) -s(index(i)) refsurf]);
            else
                disp('Refsurf must be < surface number to convert !  Please check input!');
            end
        end
    end
end

cvmacro('cvgc.seq',refsurf); 
gc = cvbuf(1,1:ims,1:6);

gc = gc(abs(s),:);

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     