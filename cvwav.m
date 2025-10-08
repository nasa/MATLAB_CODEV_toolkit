function [rmsval,xshift,yshift,strehl,weightedrms,weightedstrehl] = cvwav(nrd,nom,f,w,z)
%CVWAV gets wavefront analysis data from CodeV
%
%   function [rmsval,xshift,yshift,strehl,weightedrms,weightedstrehl] =
%               cvwav(nrd,nom,f,w,z)
%
%   NOTE: THIS FUNCTION IS REPLACED BY CVRMS
%
%   INPUTS: nrd = number of rays across diameter of pupil
%           nom = nominal "as is" condition of lens
%                   0 for best focus
%                   ~0 for "as is" (default = 1) 
%           f = field number, default = 1
%           w = wavelength number, default = 1
%           z = zoom number, default =1
%
%   OUTPUT:
%                 
%   NOTE:  tip/tilt is ALWAYS removed for this analysis
%
%   See also CVEVA, CVPMA, CVPSF, CVRSI
%

disp('This function is replaced by CVRMS and may not be supported in future releases!')
if nargin<1, nrd = 64; end
if nargin<2, nom = 1; end %default to nominal wavefront
if nargin<3, f = 1; end
if nargin<4, w = 1; end
if nargin<5, z = 1; end

if nom==1, 
    NOMstring = [' nom yes; bes no;'];
else NOMstring = [' nom no; bes yes;']; 
end
NRDstring = [' nrd ' int2str(nrd) ';'];

if nargin<3 %get all field data
    cvcmd(['buf yes; buf del b0; wav;' ...
            NOMstring ...
            NRDstring ...
           ' go; buf no;']);
    cvcmd('buf fnd b0 "FRACT"'); %jump to FRACT line in buffer
    istart = cveva(['(buf.i)'])+2;  %data starts 2 rows after "FRACT"
    [nums,numf] = cvnum;
    for i=1:numf
        xrow(i,:) = cvbuf(0,istart,2:4); %x field fraction, x field value, xshift
        yrow(i,:) = cvbuf(0,istart+1,2:7); %y field fraction, y field value, rmsval, yshift, strehl
        istart = istart+3;
    end
    lastrow = cvbuf(0,istart-1,2:3);
    xshift = xrow(:,3);
    if nom==1,                % Used for cell alignment of nominal focus output
        yshift = yrow(:,4);
        rmsval = yrow(:,3);
        strehl = yrow(:,5);
    else yshift = yrow(:,3);    % Used for cell alignment of best focus output
        rmsval = yrow(:,5);
        strehl = yrow(:,6); end
    weightedrms = lastrow(1);
    weightedstrehl = lastrow(2);
else % case where only single field/wave/zoom data is desired   
    %******* set field wave and zoom in CodeV
    [nums,numf,numw,numz] = cvnum();
    if numf>1, [fx,fy] = cvf(f); end
    if numw>1, [wl,ww] = cvw(w); end
    if numz>1, cvz(z); end
    %*******
    cvcmd(['buf yes; buf del b0; wav;' ...
            NOMstring ...
            NRDstring ...
           ' go; buf no;']);
    cvcmd('buf fnd b0 "FRACT"'); %jump to FRACT line in buffer
    istart = cveva(['(buf.i)'])+2;  %data starts 2 rows after "FRACT"
    [nums,numf] = cvnum;
    for i=1:numf
        xrow(i,:) = cvbuf(0,istart,2:4); %x field fraction, x field value, xshift
        yrow(i,:) = cvbuf(0,istart+1,2:7); %y field fraction, y field value, rmsval, yshift, strehl
        istart = istart+3;
    end
    lastrow = cvbuf(0,istart-1,2:3);
    xshift = xrow(:,3);
    if nom==1,                % Used for cell alignment of nominal focus output
        yshift = yrow(:,4);
        rmsval = yrow(:,3);
        strehl = yrow(:,5);
    else yshift = yrow(:,3);    % Used for cell alignment of best focus output
        rmsval = yrow(:,5);
        strehl = yrow(:,6); end
    weightedrms = lastrow(1);
    weightedstrehl = lastrow(2);    
    %***** reset original field and wave (ignore zoom)
    if numf>1, cvf(fx,fy); end
    if numw>1, cvw(wl,ww); end
    %*******

end
 
% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     