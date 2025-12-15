function success = cvbufput(data,buffer)
%CVBUFPUT puts array data into the CodeV buffer
%
%   function success = cvbufput(data,buffer)
%
%   INPUT:  data = 2d array of data (i.e. matrix)
%                  NOTE: if data is a scalar, then data = [data 0] to make an array
%           buffer = number of buffer, default = b1
%           
%   OUTPUT: rows, cols = size of data
%
%   See also CVBUF
%  

global CodeV

if nargin<2, buffer=1; end
if nargin<1, data = [0 0]; end

if isscalar(data), data = [data 0]; end

data = double(data);

[numrows,numcols] = size(data);
cvcmd(['buf del b' num2str(buffer)]);
% [success] = invoke(CodeV, 'ArrayToBuffer', rows, cols, buffer, data);
success = CodeV.ArrayToBuffer(numrows, numcols, buffer, data);
% success = CodeV.ARRAY_TO_BUFFER(data,buffer,0);
if success~=0, disp('Data transfer failed! '); return; end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     

return

%CODE V: Example_CVApplication.m
    % % Get the buffer data into an array. This call returns an object with the array in it.
    % cvresults = cv.BufferToArray(1, 10, 1, 10, 1);
    % disp(['BufferToArray status is ' num2str(cvresults.Item('Result'))]);
    % output = cvresults.Item('Output')
    % 
    % cvresults = cv.BUFFER_TO_ARRAY(1, 1, 10, 1, 10, 1); % Get the buffer and transpose it
    % disp(['BUFFER_TO_ARRAY status is ' num2str(cvresults.Item('Result'))]);
    % output2 = cvresults.Item('Output')
    % 
    % % Write the output from the BufferToArray call back to buffer B2
    % rc = cv.ArrayToBuffer(10,10,2, output);
    % Command(cv, 'buf lis b2');
    % 
    % % Write the generated array data to buffer B3
    % rc = cv.ArrayToBuffer(10,10,3, generatedArray);
    % Command(cv, 'buf lis b3');
    % 
    % % Write the transposed generated array data to buffer B4 
    % rc = cv.ARRAY_TO_BUFFER(generatedArray, 4, 1);
    % Command(cv, 'buf lis b4');
    % 
    % Visual Basic Syntax
    % ARRAY_TO_BUFFER(psaInput(,) As Double, bufNum As Long, transpose As Long) As Double
    % Parameters
    % psaInput	Input array.
    % bufNum	Worksheet buffer number to which the macro array is copied.
    % transpose	Whether to transpose the data in the array before copying to buffer (0 = no transpose, non-zero = transpose).
    % Return Value
    % Contains a value indicating success or failure:
	%     0	Successful copy to empty buffer.
	%     -1	Failure due to invalid inputs or some other error.
    %     n Successful copy to a buffer that was not empty. The positive integer, n, indicates the row number of the last row in the buffer that contained data before data from the array was appended.
    % 
    % Visual Basic Syntax
    % BUFFER_TO_ARRAY(bufNum, startRow, endRow, startCol, endCol, transpose) 
    % Parameters
    % bufNum	Worksheet buffer number from which numeric data will be copied.
    % startRow	Starting row number of the source area in the buffer.
    % endRow	Ending row number of the source area in the buffer
    % startCol	Starting column number of the source area in the buffer.
    % endCol	Ending column number of the source area in the buffer.
    % transpose	Whether to transpose the data in the buffer before copying to the array (0 = no transpose, non-zero = transpose).
    % Dictionary Object
    % Result dictionary item: The possible return values are:
    % 0	Data was copied successfully. For any empty or non-numeric buffer cells, corresponding destination array elements are set to 0.0.
    % 1	Buffer location specified was empty so no data was copied to array. All destination array elements are set to 0.0.
    % -1	The copy was not successful due to invalid inputs or some other error.
    % Output dictionary item: Output array
    % Example:
    % Buffer: col 1 col 2 row 1: 3 5 row 2: 4 6
    % array(1,1) = 3 array(1,2) = 5 array(2,1) = 4 array(2,2) = 6
    % Buffer and Array Methods
