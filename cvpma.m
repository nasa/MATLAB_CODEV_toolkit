function [opd,mask,data] = cvpma(TGR,fit_option,f,w,z,NRD_GRI)
%CVPMA gets CODE V exit pupil map, mask, and other PMA data
%
%   function [opd,mask,pin,Xfno,Yfno,RefRad,Xfl,Yfl] = cvpma(TGR,fit_option,f,w,z,NRD_GRI);
%
%   INPUTS: TGR: gridsize, pupil fills grid unless specific NRD_GRI option is input (default = 32)
%           fit_option: default = 0 (no fit), 1 = best fit tilt removed, 2 = best fit focus
%           f,w,z = field number, wavelength number, zoom position (all default=1)
%           NRD_GRI:  >0 = NRD input number of rays across diameter of entrance pupil
%                     default NRD = TGR input to fill array, TGR then defaults to nextpow2(TGR)
%                     <0 = GRI input, image gridsize spacing
%   OUTPUT: opd: units in meters
%           mask:  is a logical array, 0 for ray failure, 1 for pass
%           data: contains a 2 column vector of the following items:
%                 data = [wl_nm RefRad;     Wavelength (nm) and Ref Sphere Rad (cvunits)
%                         XYfield';         X and Y field
%                         XYfno';           X and Y f/#
%                         XYfl';            X and Y focal length
%                         Xentpup;          X Entrance pupil increments
%                         Yentpup;          Y Entrance pupil increments
%                         def array;        defocus value and CodeV TGR array size
%                         ];
%
%           USAGE NOTE:  GRI = Lambda_Short * F/# * NRD / TGR
%                              or, the pupil grid fills up half the
%                              transform grid for the shortest wavelength.
%                              or, AIRY disk diameter spans 4.88 grid
%                              elements at shortest wavelength
%                        "Q" = TGR/NRD = Lambda_Short * F/# / GRI
%                              =2 implies a critically sampled pupil
%                              <2 implies undersampling the pupil
%                              >2 implies oversampling the pupil
% 
%   See also cvpin, cvrac, cvhfr, cvppr, cvpsf, cvrsi, cveva, cvwav
%

buffer = 1;   % CODE V buffer used to store and transfer PMA data

if nargin<1, TGR = 32; end;
if nargin<6, NRD_GRI=TGR; end; %if no NRD_GRI input, assume the user wants output size of TGR
TGR = pow2(nextpow2(TGR));  % CodeV only accepts powers of 2 for TGR TGRs
if nargin<2, fit_option = 0; end
    if fit_option == 0, ADWstring = ' ';
    elseif fit_option == 1, ADWstring = 'ADW TIL 0 0 0 0 ;';
    elseif fit_option == 2, ADWstring = 'ADW FOC 0 0 0 0 ;';
    else ADWstring = ' ';
end
if nargin<3, f=1; end;
if nargin<4, w=1; end;
if nargin<5, z=1; end;

%******* set single field wave and zoom in CodeV
[nums,numf,numw,numz] = cvnum();
if numf>1, [fx,fy] = cvf(f); end
if numw>1, [wl,ww] = cvw(w); else [wl,ww] = cvw; end
if numz>1, cvz(z); end
%*******

if NRD_GRI>0, NRD_GRIstring = ['nrd ' int2str(NRD_GRI) ';'];
else NRD_GRIstring = ['gri ' num2str(-NRD_GRI) ';'];
end
command_string = [ ...
        'buf del b' int2str(buffer) ';' ...
        'pma;' ...
        'tgr ' int2str(TGR) ';' ...  %tgr: transform grid (16 or more, power of 2) 
        NRD_GRIstring ... %'nrd z1 ' int2str(NRD_GRI) ';' ...  %nrd: number rays across diameter NOTE: this is set to NRD to ensure an array without padding (i.e. the pupil doesn't extend to the edges)
        'lis n; ' ...  %lis n surpresses output, lis y used to be used to enable PMA output to goto buffer
        ADWstring ...
        ' wbf b' int2str(buffer) ';' ...
        'go;'];
cvcmd(command_string);

opd = cvbuf(buffer,17:(17+TGR-1),1:TGR); %opd in units of waves
wl_nm = cvbuf(buffer,6,2);
mask = opd~=-99999; % -99999 is default ray failure in opd data
opd = opd.*mask; %zero masked data
opd = opd*wl_nm*1e-9; %convert output to meters
% opd = opd.*wl(w).*1e-6.*cvunits; %This will convert opd data into lens units
% opd = opd.*wl(w).*1e-9; %This will convert opd data into meters

% Retrieve other data if requested
if nargout>2
    XYfield = cvbuf(buffer,7:8,2);
    DXDYentpup = cvbuf(buffer,12:13,2:3); %row 1 = DX, row 2 = DY
    XYfno = cvbuf(buffer,7:8,5);
    RefRad = cvbuf(buffer,9:11,2);
    XYfl = RefRad(2:3); RefRad = RefRad(1);
    def_array = cvbuf(buffer,15:16,2);
    
    data = [wl_nm RefRad;
            XYfield';
            XYfno';
            XYfl';
            DXDYentpup;
            def_array';
            ];
            
%     Xfno = cveva(['(buf.num b' int2str(buffer) ' i7 j5)']);
%     Yfno = cveva(['(buf.num b' int2str(buffer) ' i8 j5)']);
%     RefRad = cveva(['(buf.num b' int2str(buffer) ' i9 j2)']);
%     Xfl = cveva(['(buf.num b' int2str(buffer) ' i10 j2)']);
%     Yfl = cveva(['(buf.num b' int2str(buffer) ' i11 j2)']);
end

% modify data output size to input TGR request, unless specific NRD_GRI given
if nargin<6 & NRD_GRI>0 & NRD_GRI<TGR %only if NRD<TGR
    low = 1+ceil( (TGR-NRD_GRI)/2 ); %takes into account odd or even NRD
    hi = TGR - floor( (TGR-NRD_GRI)/2 );
    opd = opd(low:hi,low:hi);
    mask = mask(low:hi,low:hi);
end

% plot output if no output requested
if nargout<1
    rmsopd = norm(opd(mask)-mean(opd(mask)))./sqrt(length(find(mask))); %calculates the opd rms.
    opd(mask==0) = NaN;
    pvval = max(opd(mask))-min(opd(mask));
    figure('color',[1 1 1]), 
    if max(opd(mask))-min(opd(mask))<1e-5,
        imagesc(flipdim(opd.*1e9,1))
        title(['Exit Pupil OPD Map (Units=nm), for f,w,z']);
        xlabel(['RMS OPD = ' num2str(rmsopd*1e9) ' nm ' sprintf('\n') 'PV OPD = ' num2str(pvval*1e9) ' nm']);
    elseif max(opd(mask))-min(opd(mask))<1e-2
        imagesc(flipdim(opd.*1e6,1))
        title(['Exit Pupil OPD Map (Units=\mum)']);
        xlabel(['RMS OPD = ' num2str(rmsopd*1e6) ' \mum ' sprintf('\n') 'PV OPD = ' num2str(pvval*1e6) ' \mum']);
    else 
        imagesc(flipdim(opd.*1e3,1))
        title(['Exit Pupil OPD Map (Units=mm)']);
        xlabel(['RMS OPD = ' num2str(rmsopd*1e3) ' mm ' sprintf('\n') 'PV OPD = ' num2str(pvval*1e3) ' mm']);
    end
    colorbar, axis tight, axis square, colormap('jet')
    alpha(1*mask);
    set(gca,'XTick',[]); set(gca,'YTick',[]); set(gca,'ZTick',[])
end

%***** reset original field and wave (ignore zoom)
if numf>1, cvf(fx,fy); end
if numw>1, cvw(wl,ww); end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     