function [opl,mask] = cvopl(nrd,f,w,z,chkap)
%CVOPL uses CodeV RAYTRA engine to calculate OPL data at image surface.
%
%   function [opl,mask] = cvopl(nrd,f,w,z);
%
%   INPUTS: nrd =  number rays across diameter of entrance pupil, default = 32
%           f,w,z = number of the system positions (default is f,w,z =1)
%
%   OUTPUT: OPL = Optical path length from surface 1 to image surface in
%                   lens units (typically mm)
%           mask = Defines the pupil mask based on a ray failure map.
%           
%   See also, CVSENS, CVLOM, CVRAYGRID

if nargin<1, nrd = 32; end
if nargin<2, f=1; end 
if nargin<3, w=1; end 
if nargin<4, z=1; end
if nargin<5, chkap=1; end

cvmacro('cvopl.seq',[nrd f w z chkap]);
opl = cvbuf(1,1:nrd,1:nrd);
mask = opl~=-9999; %ray failure flag value
opl(~mask) = 0;

if nargout<1
    opd = mean(opl(mask)) - opl;
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
    colorbar, axis tight, axis square,
    set(gca,'XTick',[]); set(gca,'YTick',[]); set(gca,'ZTick',[])
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     