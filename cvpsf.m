function [psf,grid_spacing,centroid] = cvpsf(TGR,NRD_GRI,PGR,PRO,f,z)
%CVPSF  gets CodeV PSF normalized to perfect lens (i.e. strehl value)
%
%   function [psf,grid_spacing,centroid] = cvpsf(TGR,NRD_GRI,PGR,PRO,f,z);
%
%   INPUTS: TGR = Transform Grid, >16 in powers of 2, default = 128
%           NRD_GRI = Image sampling by # of rays across pupil diameter
%                   (NRD) or by actual Grid spacing at image (GRI)
%                   positive input number implies NRD, negative implies GRI
%                   default is NRD = TGR/2 (Nyquist sampling), giving
%                   image gridsize = Airy disc/4.88
%           PGR = Size of output array, default = TGR
%           PRO = (0 no,1 yes): Use propogate equations for defocused image, def = 0
%           f,z = field number, zoom position (all default=1)
%
%   OUTPUT: psf = psf data
%           grid_spacing = image grid space
%           centroid = [x y] centroid shift from chief ray coordinate at image
% 
%   See also CVEVA, CVPMA, CVRSI, CVWAV
%

CVbuf = 12;   % CODE V buffer used to store and transfer PSF data

if nargin<1, TGR = 128; end;
    TGR = pow2(nextpow2(TGR));  % CodeV only accepts powers of 2 for gridsizes
if nargin<2, NRD_GRI = TGR/2; end; % Nyquist sampling
if nargin<3, PGR = TGR-1; end;
if nargin<4, PRO = 0; end;
if nargin<5, f=1; end
if nargin<6, z=1; end

%******* set field and zoom in CodeV
[nums,numf,numw,numz] = cvnum();
if numf>1, [fx,fy] = cvf(f); end
if numz>1, cvz(z); end
%*******

TGRcmd = ['TGR ' int2str(TGR) ';'];
PGRcmd = ['PGR ' int2str(PGR) ';'];

if NRD_GRI>0, 
    NRD_GRIcmd = ['NRD ' int2str(NRD_GRI) ';']; 
    else NRD_GRIcmd = ['GRI ' num2str(-NRD_GRI) ';']; 
end;

if PRO==1, 
    PROcmd = 'PRO yes;';
    else PROcmd = ''; 
end;
command_string = ['buf del b0; buf y; buf del b' ... %del buf b0 for psf centroid data
    int2str(CVbuf) '; PSF; INT STR;' ...
    TGRcmd PGRcmd NRD_GRIcmd PROcmd ...
   'lis y; WBF b' int2str(CVbuf) '; GO; buf n;'] ; 
cvcmd(command_string);

psfgridsize = 2*ceil(PGR/2); %Output gridsize is even number
psf = cvbuf(CVbuf, 11:(11+psfgridsize-1), 1:(psfgridsize));

grid_spacing = cveva(['(buf.num b' int2str(CVbuf) ' i9 j2)']);

if nargout>2
    cvcmd(['lcl num ^i; buf fnd b0 "centroid"; ^i == (buf.i)+1; '...
            'buf put b2 i1 j1 (buf.num i^i j2); buf put b2 i1 j2 (buf.num i^i j3);']);
    centroid = cvbuf(2,1,1:2);
end
%***** reset original field (ignore zoom)
if numf>1, cvf(fx,fy); end

if nargout<1
    figure, hold on, imagesc(psf), axis equal;
    colorbar, axis tight;
    title('PSF (peak value is strehl ratio)');
    xlabel(['Image Grid Spacing = ' num2str(grid_spacing)]);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     