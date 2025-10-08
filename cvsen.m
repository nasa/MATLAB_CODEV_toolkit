function [xy,xyobj,rms,opd,mask,opd0,xydxdy] = cvsen(nrd,s,f,w,z,dx_da)
%CVSEN gets CODE V rigid body motion wavefront sensitivity data
%
%   function [xy,xyobj,rms,opd,mask,opd0,xydxdy] = cvsen(nrd,s,f,w,z,dx_da)
%
%   INPUTS:   nrd = number of rays across pupil, default = 64
%             s = surfaces that will be perturbed in rigid body
%             f w z = field wave zooom (default = 1)
%             dx_da = amount of perturbation, units in meters radians, 
%                     default = 1 um, 1 urad
%
%   OUTPUT:   sensitivities returned in m/m and m/rad = change/dx_da
%             xy = two column array of image motion of psf (based on ref sphere center)
%             xyobj = two column array of image motion of psf due to object shift
%             rms = rms value of opd from cvwav function
%             opd, mask = opd sensitivities and mask data
%             xydxdy = x y dx dy raw data in 3 dimensions
%                      row = s, col = 6 DOF, 3rd dim = x,y,dx,dy
%
%   USAGE NOTE: to get LOS data with respect to object, perform xy/xyobj
%
%   See also CVR, CVRMSWE, CVPMA, CVRBSHIFT
%

if nargin<1, nrd = 32; end;
if nargin<2, s = 1; end;
if nargin<3, f = 1; end;
if nargin<4, w = 1; end;
if nargin<5, z = 1; end;
if nargin<6, dx_da = 1e-6; end; %1 um, 1urad

if nargout<1 & length(s)>1, s = s(1); return, end %plot only one set of data

% get system data, and set single f w z in model
units = cvunits; %1 = mm, 10 = cm, 25.4 = in
[fx0,fy0] = cvf(f);
[wl0,ww0] = cvw(w); w_nm = wl0(w);
cvz(z);
ims = cvnum;

% gc0 = cvgc(1:ims); %get global data for rigid body checks

% nominal chief ray coordinates and wavefront data
r1 = cvr([1 1 ims 1 z]); x0 = r1(1); y0 = r1(2); 
[rms0,dx0,dy0] = cvrmswe(nrd,1,1,z,1);
if nargout>3 | nargout<1, [opd0,mask] = cvpma(nrd,1,1,1,z); end %best fit tilt for cvpma

% get object sensitivity data 
fldtyp = cveva('(typ fld)',1); if fldtyp == 'ANG', df = dx_da*180/pi; else df = dx_da; end
cvf(fx0(f)+df,fy0(f)); 
    % evaluate system
        r1 = cvr([1 1 ims 1 z]); x = r1(1); y = r1(2); [rms,dx,dy] = cvrmswe(nrd,1,1,z,1); %cvwav(nrd,1,1,1,z);
    % calculate sensitivities
        xyobj(1,1) = (x+dx-x0-dx0); xyobj(1,2) = (y+dy-y0-dy0);
cvf(fx0(f),fy0(f)+df); 
    % evaluate system
        r1 = cvr([1 1 ims 1 z]); x = r1(1); y = r1(2); [rms,dx,dy] = cvrmswe(nrd,1,1,z,1);%wav(nrd,1,1,1,z);
    % calculate sensitivities
        xyobj(2,1) = (x+dx-x0-dx0); xyobj(2,2) = (y+dy-y0-dy0);
cvf(fx0(f),fy0(f)); %reset field point
disp('Completed taking object data.');

% surface sensitivity data
opdsen = [];
for i=1:length(s) %loop on number of surfaces
    for j=1:6 %rigid body DOF
    % perturb surface
        cvrbshift(s(i),j,dx_da,1); %cvrbshift(s,DOFnum,shiftval,mrads)
    % evaluate system
        r1 = cvr([1 1 ims 1 z]); x(i,j) = r1(1); y(i,j) = r1(2); 
        [rms(i,j),dx(i,j),dy(i,j)] = cvrmswe(nrd,1,1,z,1);%wav(nrd,1,1,1,z);
        if nargout>3 | nargout<1, [opd(:,:,i,j),m] = cvpma(nrd,1,1,1,z); mask = mask & m; end
%     % check for independent rigid body motion (i.e. no "tail wagging" of other optical surfaces in system
%         if j<4 %only look at x y z translations
%             dgc = cvgc(1:ims)-gc0;  
%             dgc_yn = abs(dgc)>eps/2;
%             if dgc_yn(i,j)~=1, disp(['Warning! Surface ' int2str(s(i)) ' DOF ' int2str(j) ' FAILED self motion test']); dgc, end;
%             [R,C] = find(dgc_yn); surfFAIL = R(R~=s(i));
%             if length(surfFAIL)>0, disp(['Warning! Surface ' int2str(s(i)) ' DOF ' int2str(j) ' FAILED system motion test']); dgc, end; 
% %             if sum(dgc_yn(:))>1, disp(['Warning! Surface ' int2str(s(i)) ' DOF ' int2str(j) ' FAILED system motion test']); dgc, end;
%         end
    % reset surface
        cvrbshift(s(i),j,-dx_da,1);
    end
    disp(['Completed surface ' int2str(i) ' of ' int2str(length(s))]);
end

% calculate sensitivities
for i=1:length(s) %loop on number of surfaces
    for j=1:6 %rigid body DOF
        xsen(i,j) = (x(i,j)+dx(i,j)-x0-dx0);
        ysen(i,j) = (y(i,j)+dy(i,j)-y0-dy0);
        rmssen(i,j) = sqrt( abs( rms(i,j)*rms(i,j) - rms0*rms0 ) ); %root square difference
        if nargout>3 | nargout<1, opdsen(:,:,6*(i-1)+j) = mask.*(opd(:,:,i,j)-opd0); end    
    end
end

% Data for troubleshooting
% assignin('base','x',x);assignin('base','y',y);assignin('base','dx',dx);assignin('base','dy',dy);
% assignin('base','x0',x0);assignin('base','y0',y0);assignin('base','dx0',dx0);assignin('base','dy0',dy0);

xsen = xsen'; ysen = ysen'; %prepare for column formation
xy = [xsen(:) ysen(:)]*1e-3*units/dx_da; % convert to m/m and m/rad
xyobj = xyobj*1e-3*units/dx_da; % convert to m/m and m/rad
rms = rmssen*w_nm*1e-9/dx_da; %convert rms data to m/m and m/rad
if nargout>3 | nargout<1,
    opd = opdsen*w_nm*1e-9/dx_da; %convert opd data to m/m and m/rad
    opd0 = opd0*w_nm*1e-9;
end
xydxdy = cat(3,x,y,dx,dy)*1e-3*units; %convert to m/rad

% reset fps and wl
cvf(fx0,fy0);
cvw(wl0,ww0);

if nargout<1 & length(s)<2 %plot data for single case
    figure;
    for i=1:6
        data = opd(:,:,i); data(~mask) = nan;
        subplot(3,3,i), imagesc(data), axis equal, axis square, axis tight, axis off, colorbar;
        %stdval = std(data(mask)); xlabel(['rms = ' num2str(stdval)]);
    end
    subplot(3,3,1), title('dX'); subplot(3,3,2), title('dY'); subplot(3,3,3), title('dZ');
    subplot(3,3,4), title('tX'); subplot(3,3,5), title('tY'); subplot(3,3,6), title('tZ');
    subplot(3,3,7), plot(xy(:,1),'x'), hold on, plot(xy(:,2),'o'), xlim([0 7]), title('X and Y sens: {x,o}');
    subplot(3,3,8), plot(xyobj(:,1),'x'), hold on, plot(xyobj(:,2),'o'), xlim([0 3]), title('X and Y object sens: {x,o}');
    subplot(3,3,9), plot(rms(:),'o'), xlim([0 7]), title('RMS OPD');
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     