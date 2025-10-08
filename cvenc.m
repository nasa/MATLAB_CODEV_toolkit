function [radii,PER] = cvenc(PER,TGR,NRD_GRI,PGR,PRO,CEN)
%CVENC gets CODE V PSF based encircled energy for fieldpoint 1
%
%   function [radii,percent] = cvenc(PER,TGR,NRD_GRI,PGR,PRO,CEN)
%
%   INPUTS: PER = percent values to get encircled energy radius, default = 10:10:90
%           TGR = Transform Grid, >16 in powers of 2, default = 128
%       NRD_GRI = Image sampling by # of rays across pupil diameter
%                   (NRD) or by actual Grid spacing at image (GRI)
%                   positive input number implies NRD, negative implies GRI
%                   default is NRD = TGR/2 (Nyquist sampling), giving
%                   image gridsize = Airy disc/4.88
%           PGR = Size of output array, default = TGR
%           PRO = (0 no,1 yes): Use propogate equations for defocused image, def = 0
%           CEN = Centered on Chief Ray for Centroid: 0 CHI (default), 1 CEN
%
%   OUTPUT: radii = radii of data
%           PER   = percent encircled energy
% 
%   See also CVPMA, CVPSF, CVRSI, CVWAV

if nargin<1, PER = 10:10:90;    end
if nargin<2, TGR = 256;         end
    TGR = pow2(nextpow2(TGR));
if nargin<3, NRD_GRI = TGR/4;   end
if nargin<4, PGR = TGR;         end
if nargin<5, PRO = 0;           end
if nargin<6, CEN = 0;           end

numdata = length(PER);
numcalls = ceil(numdata/10); %number of calls to cvENC function required

diam = [];
for i=1:numcalls
    % build command text
    if i<numcalls
        pervals = PER( (i-1)*10+1:i*10 );
        else
        pervals = PER( (i-1)*10+1:end );
    end
        
    PERstring = [];
    for ii=1:length(pervals); PERstring = [ PERstring 'PER ' num2str(pervals(ii)) '; ' ]; end 
    if PRO == 0, PROstring = 'PRO NO '; else PROstring = 'PRO YES '; end
    if CEN == 0, CENstring = 'ENC CHI'; else CENstring = 'ENC CEN'; end

    command_string = ['buf yes; buf del b0; PSF; ' ...
                   CENstring ...
                  ';TGR ' int2str(TGR) ...
                  ';PGR ' int2str(PGR) ...
                  ';NRD ' int2str(NRD_GRI) ';' ... 
                   PERstring ...
                   PROstring ...
                  ';GO ; BUF NO;'];

    cvcmd(command_string);
    cvcmd('buf fnd b0 "PCT"'); %jump to PCT line in buffer
    istart = cveva('(buf.i)') + 2;
    iend = istart+length(pervals)-1;
    diamvals = cvbuf(0,istart:iend,2);
	diam = [diam(:); diamvals(:)];
end

radii = diam./2;

if nargout < 1,
    figure; plot(radii,PER)
    title('Encircled Energy');
    xlabel('Radius (in mm)');
    ylabel('Percent Encircled Energy');
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     