function data = cvims2exp(f,w,z)
% CVIMS2EXP converts image surface to exit pupil
% function data = cvims2exp(f,w,z)

if nargin<1, f=1; end
if nargin<2, w=1; end
if nargin<3, z=1; end

cvf(f);
cvw(w);
cvz(z);

data = cvdec(cvnum);
%********** Get Chief ray data from unperturbed system
imsRDY = cveva('(rdy si)');
cvcmd('RDY SI 0;');
img1 = cvr([1 1 cvnum 1 1]);
imgX = img1(1); imgY = img1(2); imgZ = img1(3);
imgL = img1(4); imgM = img1(5); imgN = img1(6);
refrad = cveva('(rsr z1 f1)'); %CodeV only accepts z and f qualifiers for RSR

%********** Shift ims to EXP and set proper exit pupil curvature
expX = imgX-imgL*refrad*sign(imgN);
expY = imgY-imgM*refrad*sign(imgN); 
expZ = imgZ-imgN*refrad*sign(imgN);
expA = atand(imgM/imgN); 
expB = -asind(imgL)*sign(imgN); 
expC = 0; %no clocking for now...

% [raydata,fail,fdata] = cvraytra(f,w,z);
cvshift(cvnum,[expX,expY,expZ,expA,expB,expC]);
cvcmd(['RDY SI ' num2str(refrad,20)]);
%cvcmd(['THI SI ' num2str(-refrad,20)]);


data = [data imsRDY];
cvcmd('DEL APE SI');

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     