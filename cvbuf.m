function [buf_num,buf_cell,buftyp] = cvbuf(buf_cv,rows,cols)
% CVBUF gets buffer contents from CODE V
%
%   function [buf_num,buf_cell] = cvbuf(buf_cv,rows,cols,fail_value)
%
%   INPUT:  buf_cv = CODE V buffer number, default = 0 (command window)
%                    if input < 0, then no call to CodeV will be made
%           rows = rows desired, MUST BE CONTINUOUS (e.g. 1:4)
%           cols = cols desired, MUST BE CONTINUOUS (e.g. 3:7)
% 
%   OUTPUT: 
%           buf_num  = array of buffer data, "STR" and "UNKNOWN" cells are set to NaN
%           buf_cell = cell array of buffer data, taken from exported file (delimited by \t)
% 
%   See also CVBUFLIS, CVBUFSTR
%
%   USAGE: single row and column request returns NUM and STR using CVEVA
%          multi-row and column requests only NUM data (for now)
%                full buffer request returns NUM and STR in buf_cell output

global CodeV
if nargin<1, buf_cv=0; end
if nargin<2, rows=[]; end
if nargin<3, cols=[]; end

% if isempty(filename), filename = [cvpath '\temp.txt']; end

if isempty(rows) || isempty(cols), [rownum,colnum] = cvbufsize(buf_cv); end
if isempty(rows), rows = 1:rownum; end
if isempty(cols), cols = 1:colnum; end

% if buf_cv == 0 %move buffer 0 to empty buffer to extract data
%     buf_empty = 1;
%     while ~cvdb(['buf.emp b' int2str(buf_empty)])
%         buf_empty = buf_empty + 1;
%     end
%     buf_cv = buf_empty; %disp(buf_cv)
%     cvcmd(['buf mov b' int2str(buf_cv) ' i1 j1; buf cop b0 ia ja']);
% end

if numel(rows)==1 && numel(cols)==1 %single value requested, use cveva
    buftyp = cveva(['(buf.typ b' int2str(buf_cv) ' i' int2str(rows) ' j' int2str(cols) ')'],1);
    if contains( buftyp , 'NUM' )
        cmd = ['(buf.num b' int2str(buf_cv) ' i' int2str(rows) ' j' int2str(cols) ')'];
        buf_num = cveva(cmd,0);
    elseif contains( buftyp , 'STR' )
        cmd = ['(buf.str b' int2str(buf_cv) ' i' int2str(rows) ' j' int2str(cols) ')'];
        buf_num = cveva(cmd,1);
    else
        buf_num = buftyp;
    end
    buf_cell = {buf_num};
elseif nargin>1 %multiple values requested with specific rows and cols, use BufferToArray
    if nargout<2 %fast BufferToArray transfer for NUM only output
        cvresults = CodeV.BufferToArray(rows(1), rows(end), cols(1), cols(end), buf_cv);
        success = cvresults.Items{1};
        if success~=0
            disp('Data transfer failed!  Possibly due to non-numerical data in CV buffer.');
            buf_num = 0;
            return;
        end    
        buf_num = cell2mat( cvresults.Item('Output') );
    else %slower loop EVA transfer includes STR and NUM output
        [buf_cell,buftyp] = cvbufstr( buf_cv, rows, cols);
        buf_num = zeros(size(buftyp));
        buf_num( buftyp == 0 ) = buf_cell{buftyp == 0};
        buf_cell( buftyp == -1 ) = {[]};
    end
else %extract full buffer using export and readcell
    buftyp = cvbuftyp( buf_cv , rows, cols );
    buf_cv_str = int2str(buf_cv);
    cv_cd = cvdb('cd',1);   % Get CODE V working directory, needed when MATLAB CD ~= CODE V CD
    filename = fullfile( cv_cd, ['buf' buf_cv_str '.txt'] );
    cmd = ['buf exp b' buf_cv_str ' "' filename '"'];
    cvcmd(cmd);
    pause(0.2); %pause to write file
    buf_cell = readcell(filename,'EmptyLineRule','read','Delimiter','\t');
    delete(filename);
    buf_cell( buftyp == -1 ) = {[]};
    buf_num = zeros(size(buftyp));
    buf_num ( buftyp == 0 ) = cell2mat( buf_cell( buftyp == 0 ) );
end

if nargout<1
    disp(buf_num);
    buf_num = [];
end


% Authors: Joseph M. Howard, 
% Revision Date: 2020.02.28


% CODE V HELP
%     [double, SafeArray Pointer(double)]	BUFFER_TO_ARRAY	(handle, int32, SafeArray Pointer(double), int32, int32, int32, int32, int32)
% Visual Basic Syntax
% BUFFER_TO_ARRAY(buf_num As Integer, array() As Double, start_row As Integer, end_row As Integer, start_col As Integer, end_col As Integer, transpose As Integer) As Double
%   Parameters
% buf_num array
% start_row end_row start_col end_col transpose
% Return Value
% Worksheet buffer number from which numeric data will be copied.
% Array of doubles to be filled with the worksheet buffer data. This array must be (end_row - start_row + 1) x (end_col - start_col + 1). Starting row number of the source area in the buffer.
% Ending row number of the source area in the buffer
% Starting column number of the source area in the buffer. Ending column number of the source area in the buffer.
% Whether to transpose the data in the buffer before copying to the array (0 = no transpose, non-zero = transpose).
%  The possible return values are:
%  0  Data was copied successfully. For any empty or non-numeric buffer cells, corresponding destination array elements are set to 0.0.
%  1  Buffer location specified was empty so no data was copied to array. All destination array elements are set to 0.0.
% -1  The copy was not successful due to invalid inputs or some other error.
   
