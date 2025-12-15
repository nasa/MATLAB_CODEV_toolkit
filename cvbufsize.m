function [rows,cols] = cvbufsize(bufnum)
%CVBUFSIZE gets buffer size from CODE V
%
%   function [rows,cols] = cvbufsize(bufnum)
%
%   INPUT:  buffer = number of buffer, default = b1
%   OUTPUT: [rows,cols] = number of rows and cols of buffer
% 
%   See also CVBUF CVDB CVEVA
 
if nargin<1, bufnum = 1; end

buf_str = int2str(bufnum);
rows = cvdb(['BUF.LST B' buf_str]);
cols = cvdb(['BUF.MXJ B' buf_str]);

if nargout<1
    disp(['Size of Buffer ' int2str(bufnum)]);
    disp(['Rows, Cols: ' int2str([rows cols])]);
end

% Authors: Joseph M. Howard
% Revision Date: 2020.02.24