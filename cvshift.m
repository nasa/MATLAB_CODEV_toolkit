function cvshift(arg1, arg2, arg3, arg4)
%CVSHIFT  "shifts" the decenters relative to the current value in CODE V 
%
%   function cvshift(arg1, arg2, arg3, arg4)
%
%   INPUTS: 
%       For a single degree of freedom, any number of surfaces:
%           arg1 = surfnum: surface number to be changed
%           arg2 = DOFnum: (1-9) xde, yde, zde, ade, bde, cde, xod, yod,
%                           zod respectively
%           arg3 = shiftval, UNITS: lens units (e.g. mm), and degrees
%           arg4 = mrads: i.e. interpret input as meters radians
%                     0 = input data is assumed to be lens units
%                     1 = input data is assumed as meters radians
% 
%       For two or more degrees of freedom, any number of surfaces:
%       CASE 1: USING ONLY ONE ARGUMENT: arg1
%           arg1 = [surfnum shiftval]
%           First column contains the surface numbers to be perturbed
%           Columns 2-10 contain the shiftval, in lens units. You only need
%           to specify out as many columns are there are data. I.E.: For
%           only x- and y- shifts, you only need 3 columns: [surface; XDE;
%           YDE]
%
%       CASE 2: USING TWO ARGUMENTS: arg1, arg2
%           arg1 = surfnum: surface number to be changed (vector)
%           arg2 = shiftvals: matrix of values to be shifted. As with Case
%           1, you only need to use the number of matrix columns you need.
%           Columns are in this order: 1=XDE, 2=YDE, 3=ZDE, 4=ADE, 5=BDE,
%           6=CDE, 7=XOD, 8=YOD, 9=ZOD
%         OR
%           arg1 = [surfnum shiftvals]: surface number to be changed (vector)
%           arg2 = mrads
%
%       CASE 3: USING THREE ARGUMENTS: arg1, arg2, arg3
%           arg1 = surfnum: surface number to be changed (vector)
%           arg2 = shiftvals
%           arg3 = mrads
%           shiftval columns: [XDE YDE ZDE ADE BDE CDE XOD YOD ZOD]
%               NOTE: single data for arg1 arg2 and arg3 defaults to single DOF
%               case above
%   NOTE: These are deltas from existing values within CodeV
%
%   See also CVRBSHIFT
%

if nargin < 1, disp('Need input data!'); help cvshift, return, end;

if nargin == 1, %surface and shiftval data in arg1
    sw = 'allDOF';
    [numsurfs,c] = size(arg1);
    s = arg1(:,1);
    if c < 2, disp('Need shift data'), return, end;
    shiftval = arg1(:,2:end);
    mrads = 0; 
elseif nargin == 2,
    sw = 'allDOF';
    [numsurfs,c] = size(arg1);
    s = arg1(:,1);
    if c>2 %surface and shift data in arg1, mrads in arg2
        shiftval = arg1(:,2:end);
        mrads = arg2;
    else %surface data in arg1, shiftval data in arg2
        shiftval = arg2;
        mrads = 0;
    end
elseif nargin == 3,
    s = arg1;
    numsurfs = length(s);
    if ( length(arg2)>1 && length(arg3)==1 ) %shiftval in arg2, mrads in arg3
        sw = 'allDOF';
        shiftval = arg2;
        mrads = arg3;
    else  
        sw = 'singleDOF';
        DOFnum = arg2; 
        shiftval = arg3;
        mrads = 0; %default to lens units
    end
elseif nargin == 4,
    sw = 'singleDOF';
    s = arg1;
    numsurfs = length(s);
    DOFnum = arg2; % Replaces cvrbshift
    shiftval = arg3;
    mrads = arg4;
end

switch sw
%**************************************************************************
case 'allDOF'

[r,c] = size(shiftval);
shiftval = [shiftval zeros(r,9-c)]; %pad shiftval to full 9 columns

if mrads
    shiftval(:,[1:3 7:9]) = shiftval(:,[1:3 7:9])*1000/cvunits;
    shiftval(:,[4:6]) = shiftval(:,[4:6])*180/pi;
end
     
XDE = shiftval(:,1); YDE = shiftval(:,2); ZDE = shiftval(:,3);
ADE = shiftval(:,4); BDE = shiftval(:,5); CDE = shiftval(:,6);
XOD = shiftval(:,7); YOD = shiftval(:,8); ZOD = shiftval(:,9);

for ii=1:numsurfs
    sn = num2str(s(ii));
    if (XOD(ii)==0) & (YOD(ii)==0) & (ZOD(ii)==0),
        ODstr = ' ';
    else ODstr = ['xod s' sn ' (xod s' sn ')+' num2str(XOD(ii),20) ...
         '; yod s' sn ' (yod s' sn ')+' num2str(YOD(ii),20) ...
         '; zod s' sn ' (zod s' sn ')+' num2str(ZOD(ii),20) ];
    end
    cmd = ['xde s' sn ' (xde s' sn ')+' num2str(XDE(ii),20) ...
         '; yde s' sn ' (yde s' sn ')+' num2str(YDE(ii),20) ...
         '; zde s' sn ' (zde s' sn ')+' num2str(ZDE(ii),20) ...
         '; ade s' sn ' (ade s' sn ')+' num2str(ADE(ii),20) ...
         '; bde s' sn ' (bde s' sn ')+' num2str(BDE(ii),20) ...
         '; cde s' sn ' (cde s' sn ')+' num2str(CDE(ii),20) ...
         ';' ODstr ];
    cvcmd(cmd);
end

%**************************************************************************
case 'singleDOF'

    if mrads
        if DOFnum<4 | DOFnum>6, shiftval = shiftval*1000/cvunits;
        else shiftval = shiftval*180/pi;
        end
    end

for jj=1:numsurfs
    for ii=1:length(DOFnum)
        if     DOFnum(ii)==1, DOF = 'xde';
        elseif DOFnum(ii)==2, DOF = 'yde';
        elseif DOFnum(ii)==3, DOF = 'zde';
        elseif DOFnum(ii)==4, DOF = 'ade';
        elseif DOFnum(ii)==5, DOF = 'bde';
        elseif DOFnum(ii)==6, DOF = 'cde';
        elseif DOFnum(ii)==7, DOF = 'xod';
        elseif DOFnum(ii)==8, DOF = 'yod';
        elseif DOFnum(ii)==9, DOF = 'zod';
        else disp('Did not shift surface(s)!'); return;
        end;
        cmd = [ DOF ' s' num2str(s(jj)) ' ((' DOF ' s' ...
            num2str(s(jj)) ')+' num2str(shiftval(jj),20) ');' ];
        cvcmd(cmd);
    end
end
end   

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     