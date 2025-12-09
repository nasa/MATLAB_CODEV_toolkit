function [int,pha,gridspace,wl] = cvbpr(s,so,tgr,nrd_gri,pgr,f,z)
%CVBPR  gets CodeV beam propagation data
%   THIS FUNCTION IS STILL UNDER CONSTRUCTION!!!!
%   function [int,pha,gridspace,wl] = cvbpr(s,so,tgr,nrd_gri,pgr,f,z)
%
%   INPUTS: s = surfaces which to perform beam propagation between
%           so = output surface which to return amplitude and phase data
%           tgr = Transform Grid, >16 in powers of 2, default = 128
%           nrd_gri = Image sampling by # of rays across pupil diameter
%                   (NRD) or by actual Grid spacing at image (GRI)
%                   positive input number implies NRD, negative implies GRI
%                   default is NRD = tgr/2 (Nyquist sampling), giving
%                   image gridsize = Airy disc/4.88
%           pgr = Size of output array, default = tgr (currently not coded)
%           f,z = field number, zoom position (all default=1)
%
%   OUTPUT: int = intensity data at surface (2d stack)
%           pha = phase data at surface (units in nm)
%           gridspace = grid space of output data
%           wl = wavelength of returned data
% 
%   See also CVPSF, CVPMA
%

[nums,numf,numw,numz] = cvnum();

if nargin<1, s = nums-1:nums; end %default to final thickness
if nargin<2, so = nums; end
if nargin<3, tgr = 128; end;
    tgr = pow2(nextpow2(tgr));  % CodeV only accepts powers of 2 for gridsizes
if nargin<4, nrd_gri = tgr/2; end; % Nyquist sampling
if nargin<5, pgr = tgr-1; end;
if nargin<6, f=1; end
if nargin<7, z=1; end

%******* set field and zoom in CodeV
if numf>1, [fx,fy] = cvf(f); end
if numz>1, cvz(z); end
%*******

tgrcmd = ['tgr ' int2str(tgr) ';'];
pgrcmd = ['pgr s' int2str(s(1)) '..' int2str(s(end)) ' ' int2str(pgr) ';'];
if nrd_gri>0, 
    nrd_gricmd = ['NRD ' int2str(nrd_gri) ';']; 
    else nrd_gricmd = ['GRI ' num2str(-nrd_gri) ';']; 
end;
dbpcmd = ['dbp s' int2str(s(1)) '..' int2str(s(end)) ';'];
wbfcmd = ['wbf int sur s' int2str(so(1)) ' b1; wbf pha sur s' int2str(so(1)) ' b2;'];  %write data to buffer b1 b2

% buf del b1; bpr; nrd 64; dbp s33..i; wbf cmp sur si b1; go
command_string = ['buf del b1; buf del b2; bpr; ' ...
    dbpcmd tgrcmd nrd_gricmd pgrcmd wbfcmd...
   'go;'] ;
cvcmd(command_string);

gridsize = cvbuf(1,20,2); %Output gridsize is sometimes changes from input
int = cvbuf(1, 21:(21+gridsize-1), 1:gridsize); 
gridsize = cvbuf(2,20,2); %Output gridsize is sometimes changes from input
pha = cvbuf(2, 21:(21+gridsize-1), 1:gridsize); 
wl = cvbuf(2,3,6);
pha = pha*wl*1e-9/100; %convert phase data to nm

gridspace = cvbuf(1,19,2);
% gridspace = cvbuf(2,19,2);


%***** reset original field (ignore zoom)
if numf>1, cvf(fx,fy); end

if nargout<1
    figure, hold on, 
        subplot(1,2,1), imagesc(int), axis equal; colorbar, axis tight;
        title('INT');
        xlabel(['Grid Spacing = ' num2str(gridspace)]);
        subplot(1,2,2), imagesc(pha), axis equal; colorbar, axis tight;
        title('PHA');
        xlabel(['Units = nm']);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     