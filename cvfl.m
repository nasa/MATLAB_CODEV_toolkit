function    [fl,xyfl,dxdy] = cvfl(f,w,z,dx_da)
%CVFL evaluates the focal length
%
%   function    [fl,xyfl,dxy] = cvfl(f,w,z,dx_da);
%
%   INPUT:  f,w,z: field wave and zoom for calculation
%           dx_da: field perturbation for calcuation, default = 1nm or 1nr
%
%   OUTPUT: fl = mean focal length in x and y (units = meters)
%           xyfl = actual x and y focal lengths (units = meters)
%           dxdy = 2x2 matrix of differential ray data
%                   row 1 = [x-x0 y-y0] for dfx
%                   row 2 = [x-x0 y-y0] for dfy
%
%   See also: CVR
%

disp('Function under construction... please use another.'); return

global CodeV

if nargin<1, f=1; end
if nargin<2, w=1; end
if nargin<3, z=1; end
if nargin<4, dx_da = 1e-6/cvunits; end%1 nr or nm motion

fldtyp = cveva('(typ fld)',1); 
if fldtyp == 'ANG', 
    df = dx_da*180/pi; 
else 
    df = dx_da; 
end
units = cvunits;

[fx,fy] = cvf(f); %get field data, and set single field point

fx0 = fx(f); fy0 = fy(f);

input(3)=tand(fx0);
input(4)=tand(fy0);
input(1:2)=[-x,-y];
[success,ray0]=invoke(CodeV,'RAYTRA',z,w,0,input',zeros(8,1));

input(3)=tand(fx0+df);
input(4)=tand(fy0);
input(1:2)=[-x,-y];
[success,rayX]=invoke(CodeV,'RAYTRA',z,w,0,input',zeros(8,1));

input(3)=tand(fx0);
input(4)=tand(fy0+df);
input(1:2)=[-x,-y];
[success,rayY]=invoke(CodeV,'RAYTRA',z,w,0,input',zeros(8,1));

cvf(fx,fy); %reset field points

Xfl = ( ( rayX(1)-ray0(1) )^2 + ( rayX(2)-ray0(2) )^2 )^(1/2);
Yfl = ( ( rayY(1)-ray0(1) )^2 + ( rayY(2)-ray0(2) )^2 )^(1/2);
Xfl = Xfl*1e-3*units/dx_da;
Yfl = Yfl*1e-3*units/dx_da;

fl = mean([Xfl,Yfl]);
xyfl = [Xfl; Yfl];

dxdy = [ rayX(1)-ray0(1) rayX(2)-ray0(2); ...
         rayY(1)-ray0(1) rayY(2)-ray0(2);]*1e-3*units/dx_da; 

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     