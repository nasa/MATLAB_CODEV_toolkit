function [X,Y,Z,L,M,N,AOI,AOR,LEN] = cvrsi(f,r,s,w,z,global_coord)
%CVRSI gets CODE V ray trace data from RSI command
%
%   function [X,Y,Z,L,M,N,AOI,AOR,LEN] 
%                  = cvRSI(r,r,s,w,z,global_coord);
%
%   INPUTS:  f = field, number or relative [fx fy]
%            r = pupil, ray number or relative pupil
%            s = surface numbers for ray info, must be continuous
%                 default = all
%            w = wavelength number, default = 1
%            z = zoom position, default = 1
%            global_coord: 0=no (default), or global reference surface 
%
%   OUTPUT: 
%
%   See also CVEVA, CVPMA, CVPSF, CVWAV
%

%            relpup_relfield = [x_pupil y_pupil xfield y_field] coordinates for raytrace

if nargin<1, f = 1; end
if nargin<2, r = 1; end
if nargin<3, s = 1:cvnum; end
if nargin<4, w = 1; end
if nargin<5, z = 1; end
if nargin<6, global_coord = 0; end

if global_coord == 0
    glo_str = [' GLO N; '];
else 
    glo_str = [' GLO S' num2str(global_coord) ' ; '];
end

format = 'ROF ''E21.15'' X Y Z L M N AOI AOR LEN;';
cmd_s = [' s' int2str(s(1)) '..' int2str(s(end))];
cmd_w = [' w' int2str(w(1))];
cmd_z = [' z' int2str(z(1))];
if length(r)>1
    if length(f)>1
        cmd_fr = [num2str(r) num2str(f)];
    else
        cmd_fr = [' f' int2str(f) ' ' num2str(r)];
    end
else
   cmd_fr = [' f' int2str(f) ' r' int2str(r)]; 
end

cmd = [glo_str format 'buf del b0; buf;' ...
                  'rsi ' cmd_s cmd_w cmd_z cmd_fr ';' ...
                  'buf n;'];
cvcmd(cmd);
cvmacro('cvrsi.seq'); %parse data for extraction from buffer

cvcmd('buf fnd b0 "Position"'); %jump to Position line in buffer
istart = str2num(cveva('(buf.i)',1));
if global_coord == 0,
    istart = istart+2;
else istart = istart+3;
end

cvcmd('buf fnd b0 "IMG"'); %jump to IMG line in buffer
iend  = str2num(cveva('(buf.i)',1)); 

[success,output] = cvbuf(0,istart:iend,2:10);

X = output(:,1);   Y = output(:,2);   Z = output(:,3);
L = output(:,4);   M = output(:,5);   N = output(:,6);
AOI = output(:,7); AOR = output(:,8); LEN = output(:,9);

if nargout<1
   X = output; 
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     