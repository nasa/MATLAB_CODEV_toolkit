function [string,buftyp] = cvbufstr(buffer,rows,cols)
%
%CVBUFSTR gets a string data from the applicable CodeV buffer
%
%   function number = cvbufstr(buffer,rows,cols)
%
%   INPUT:  buffer = number of buffer, default = b0
%           rows, cols = index of buffer
%   OUTPUT: string, if single buffer cell or row is requested
%           cell array, if multiple buffers are requested
%
%   See also CVBUF, CVEVA, CVBUFTYP
%  

if nargin<1, buffer=0; end
if nargin<2, rows = 1:cveva(['(buf.lst b' num2str(buffer) ')'],0); end %returns number of rows
if nargin<3, cols = 1:cveva(['(buf.mxj b' num2str(buffer) ')'],0); end %entire buffer row is retrieved

%     if length(rows)<2
%         string = cveva(['(buf.str b' int2str(buffer) ' i' int2str(rows) ')'],1);
%     else
%         for i=1:length(rows)
%             string{i} = cveva(['(buf.str b' int2str(buffer) ' i' int2str(rows(i)) ')'],1);
%         end
%     end
% else

numr = length(rows);
numc = length(cols);
buftyp = zeros(numr,numc);
for i=1:numr
    for j=1:numc
        bij_str = ['b' int2str(buffer) ' i' int2str(rows(i)) ' j' int2str(cols(j)) ')'];
        buftypstr = cveva(['(buf.typ ' bij_str],1);
        if contains( buftypstr , 'STR' )
            buftyp(i,j) = 1;
            string{i,j} = cveva(['(buf.str ' bij_str],1);
        elseif contains( buftypstr , 'NUM' )
            buftyp(i,j) = 0;
            string{i,j} = cveva(['(buf.num ' bij_str],1);
        else
            buftyp(i,j) = -1;
            string{i,j} = '';
        end
    end
end
  
%
%   Copyright 2003-2005: NASA/GSFC, Optics Branch
%   CodeV verion 9.60
%
%   Author: Joseph M. Howard
%   Revisions: 
%   Revision: 1.0.0
%   Date: 2005/07/01
%  