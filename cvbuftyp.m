function buf_typ = cvbuftyp(buf_cv,rows,cols)
%CVBUFTYP gets buffer cell types from CODE V via CVEVA or buffer file export
%
%   function buf_typ = cvbuftyp(buf_cv,rows,cols)
%
%   INPUT:  buffer = number of buffer, default = b1
%   OUTPUT: array with the following data
%           -1 = UNKNOWN type cell
%            0 = NUM type cell
%            1 = STR type cell
%   See also CVBUF
 
if nargin<1, buf_cv = 1; end
if nargin<2, rows = []; end
if nargin<3, cols = []; end

[rownum,colnum] = cvbufsize(buf_cv);
if isempty(rows), rows = 1:rownum; end
if isempty(cols), cols = 1:colnum; end

%find an empty buffer to put typ data into
bufget = 11;
while ~cveva(['(buf.emp b' int2str(bufget) ')'])
    bufget = bufget + 1;
end

cvmacro('cvbuftyp.seq',[ buf_cv bufget ]);

cvcmd(['buf exp b' int2str(bufget) ' "buftyp.txt"; buf del b' int2str(bufget)]);
cv_cd = cvdb('cd',1);  % Get CODE V working directory, needed when MATLAB CD ~= CODE V CD
buf_file = fullfile(cv_cd,'buftyp.txt');
buf_typ = importdata( buf_file );
delete( buf_file );
% buf_typ = cvbuf(bufget,1:rownum,1:colnum); %INFINITE LOOP: cvbuf calls cvbuftyp

buf_typ = buf_typ(rows,cols);

% Authors: Joseph M. Howard
% Revision Date: 2020.02.24