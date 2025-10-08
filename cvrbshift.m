function gc = cvrbshift(s,DOFnum,shiftval,mrads)
% CVRBSHIFT "shifts" the decenters and tilts for single surface 'surfnum', 
%            relative to the current value already in CodeV
%
%   function success = cvrbshift(s,DOFnum,shiftval,mrads)
%
%   INPUTS: s: surface number to be changed
%           DOFnum: (1-9) xde, yde, zde, ade, bde, cde, xod, yod, zod respectively
%           shiftval UNITS: lens units (e.g. mm), and degrees
%           mrads: interpret input as meters and radians, default = 0 (no)
% 
%   See also CVSHIFT
%

if nargin<1, s = 1; end;
if nargin<2, DOFnum = 1; end;
if nargin<3, shiftval = 1e-6; end %1nm
if nargin<4, mrads = 0; end %1nr for angles

if mrads
    if DOFnum<4 | DOFnum>6, shiftval = shiftval*1000/cvunits;
    else shiftval = shiftval*180/pi;
    end
end

for jj=1:length(s)
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
        else return;
        end;
        cmd = [ DOF ' s' int2str(s(jj)) ' ((' DOF ' s' int2str(s(jj)) ')+' num2str(shiftval(jj),20) ');' ];
        cvcmd(cmd);
    end
end

if nargout>0, gc = cvgc(s); end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     