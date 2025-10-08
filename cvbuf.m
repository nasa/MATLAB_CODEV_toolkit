function data = cvbuf(buffer,rows,cols)
%CVBUF gets buffer contents from CODE V
%
%   function text = cvbuf(buffer,rows,cols)
%
%   INPUT:  buffer = number of buffer, default = b0
%           rows = rows desired, MUST BE CONTINUOUS (e.g. 1:4)
%           cols = cols desired, MUST BE CONTINUOUS (e.g. 3:7)
%   OUTPUT: 
%           Case 1:  nargin == 1, text array of buffer
%           Case 2:  nargin == 2, text array of buffer lines
%           Case 3:  nargin == 3, numerical array output (fails if text encountered)
% 
%   See also CVIN 

global CodeV
if nargin<1, buffer=0; end

if nargin<2 %buffer only input, export buffer to a file and read into matlab
    filename = [cvpath '\temp.txt'];
    cvcmd(['buf exp b' int2str(buffer) ' ' filename]);
    pause(0.5); %need pause to finish writing file
    fid = fopen(filename);
    linenumber = 1;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, 
        else text{linenumber}=tline; linenumber = linenumber+1;
        end
    end
    fclose(fid);
    delete(filename);
    data = strvcat(text');
    if nargout<1
        disp(data);
        data = [];
    end
elseif nargin<3 % buffer and row input, reads the entire row into matlab
    for i=1:length(rows)
        data{i} = cveva(['(buf.txt b' int2str(buffer) ' i' int2str(rows(i)) ')'],1);
    end
    data = strvcat(data');
else %buffer rows and cols given
    single = 0;
    if length(rows)==1 & length(cols)==1, single = 1; cols = [cols cols+1]; end  %sets two data elements required for transfer
    [success,data] = invoke(CodeV, 'BufferToArray', rows(1), rows(end), cols(1), cols(end), buffer, zeros(length(rows),length(cols)));
    if success~=0, disp('Data transfer failed!  Possibly due to non-numerical data in CV buffer.'); return; end
    if single==1, data = data(1); end %remove extra data if only one desired
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     