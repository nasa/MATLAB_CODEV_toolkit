function [opd,mask,X,Y,Z,L,M,N] = cvexp(nrd,f,w,z)
%CVEXP uses CodeV RAYTRA engine to calculate ray data at exit pupil
%
%   function [OPL,mask,X,Y,Z,L,M,N] = cvexp(nrd,f,w,z);
%
%   INPUTS: nrd =  pupil rayset density, default = 32
%           f,w,z = number of the system positions (default is f,w,z =1)
%
%   OUTPUT: OPD = Optical path difference in lens units (typically mm)
%           mask = Defines the pupil mask based on a ray failure map.
%           X,Y,Z,L,M,N = grid of data for given parameter.
%           
%   See also, CVSENS, CVLOM, CVRAYGRID, CVOPL

if nargin<1, nrd = 32; end;
if nargin<2, f=1; end; 
if nargin<3, w=1; end; 
if nargin<4, z=1; end

ims2expdata = cvims2exp(f,w,z); %create dummy exit pupil

if nargout<3,
    [opl,mask] = cvopl(nrd,f,w,z); %faster routine for opl data only
else
    [data,mask] = cvrgrid(nrd,f,cvims,w,z);
    opl = data(:,:,7);
    X = data(:,:,1); Y = data(:,:,2); Z = data(:,:,3);
    L = data(:,:,4); M = data(:,:,5); N = data(:,:,6);
end
opd = mean(opl(mask))-opl; opd(mask==0) = 0;
    
cvexp2ims(ims2expdata); %return dummy exit pupil to image surface

if nargout<1 
    rmsopd = norm(opd(mask)-mean(opd(mask)))./sqrt(length(find(mask))); %calculates the opd rms.
    opd(mask==0) = NaN;
    pvval = max(opd(mask))-min(opd(mask));
    figure,
    if max(opd(mask))-min(opd(mask))<1e-5,
        imagesc(flipdim(opd.*1e6,1))
        title(['Exit Pupil OPD Map (Units=nm)']);
        xlabel(['RMS OPD = ' num2str(rmsopd*1e6) ' nm ' sprintf('\n') 'PV OPD = ' num2str(pvval*1e6) ' nm']);
    elseif max(opd(mask))-min(opd(mask))<1e-2
        imagesc(flipdim(opd.*1e3,1))
        title(['Exit Pupil OPD Map (Units=\mum)']);
        xlabel(['RMS OPD = ' num2str(rmsopd*1e3) ' \mum ' sprintf('\n') 'PV OPD = ' num2str(pvval*1e3) ' \mum']);
    else 
        imagesc(flipdim(opd,1))
        title(['Exit Pupil OPD Map (Units=mm)']);
        xlabel(['RMS OPD = ' num2str(rmsopd) ' mm ' sprintf('\n') 'PV OPD = ' num2str(pvval) ' mm']);
    end
    colorbar, axis tight, axis square
    set(gca,'XTick',[]); set(gca,'YTick',[]); set(gca,'ZTick',[])
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     