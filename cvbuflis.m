function [buf_lis,buf_lis_cell] = cvbuflis(buf_cv,rows,NOL)
% CVBUFLIS gets listed buffer contents from CODE V, displayed precision
%
%   function [buf_lis,buf_lis_cell] = cvbuflis(buf_cv,rows,NOL)
%
%   INPUT:  buf_cv = CODE V buffer number, default = 0 (command window)
%                    if input < 0, then no call to CodeV will be made
%           rows = rows desired, MUST BE CONTINUOUS (e.g. 1:4)
%           NOL = "true" to request no lines in buf_lis output
%   OUTPUT: 
%           buf_lis  = char list of contents in buffer
%           buf_lis_cell = cell vector of buf_lis 
% 
%   See also CVBUF

if nargin<1, buf_cv=0; end
if nargin<2, rows = []; end
if nargin<3, NOL = false; end

if NOL, NOLstr = ' NOL'; else NOLstr = []; end
buf_lis = cvcmd(['buf lis b' int2str(buf_cv) NOLstr]);
buf_lis_cell = splitlines(buf_lis); %reformat to char array
buf_lis = char( buf_lis_cell );

if ~isempty(rows)
    buf_lis_cell = buf_lis_cell(rows);
    buf_lis = buf_lis(rows,:);
end

if nargout<1
    disp(buf_lis);
end


% Authors: Joseph M. Howard, 
% Revision Date: 2020.02.28

