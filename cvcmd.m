function [output,cmd,output_cell] = cvcmd(cmd)
%CVCMD  sends command to CODE V command line
%
%   function  [output,cmd] = cvcmd(cmd)
%
%   cmd = CodeV command line script to be executed, can be cell vector 
%         NOTE: multiple lines are written to a temp file and executed
%         using cvin.
%   Output is text returned from CodeV command line
%
%   See also CVON, CVOFF, CVDB, CVIN, CVOPEN, CVSAVE
%

global CodeV

output = [];
output_cell = {};
if isempty(cmd), return; end

size_cmd = size(cmd);

if size_cmd(1)<2  %single row input
    if size_cmd(2)<257
%         output = invoke(CodeV,'Command',cmd); %pre v11.3
        output = CodeV.Command(cmd);
    else
        if nargout>0, disp('Writing temporary command file.'); end
        fid = fopen('cmd.seq','wt','native','US-ASCII');
        fprintf(fid,'%s',cmd');
        fclose(fid);
        [~,output] = cvin('cmd.seq'); %load temp .seq file into CODE V
        delete cmd.seq;
    end
else  %array input, write a file and read into CODE V
    
    fid = fopen('cmd.seq','wt','native','US-ASCII');
    
    if ischar(cmd) %char array
        for i=1:size(cmd,1)
            fprintf(fid,'%s\n',cmd(i,:));
        end
    elseif isstring(cmd) %string array
        for i=1:length(cmd(:))
            fprintf(fid,'%s\n',cmd(i,:));
        end        
    elseif iscell(cmd) %string array
        for i=1:length(cmd)
            fprintf(fid,'%s\n',cmd{i});
        end        
    end
    
    fclose(fid);
    [~,output] = cvin('cmd.seq'); %load temp .seq file into CODE V
    delete cmd.seq;
    
end


% remove all "carriage return" characters (=13), and only keep "new line"
% characters (=10)
output = output( output~=13 );

% create output in cell format, which is useful for export into excel using
% xlswrite
newlines = find( output == 10 );
for i=1:length(newlines)-1
    output_cell{i,1} = output( newlines(i)+1:newlines(i+1)-1 );
end


%ARCHIVE CODE
%     %required by some unknown bug
%     fid = fopen('cmd.seq','w');
%     for i=1:numrows
%         cmdline = cmd(i,:);
%         
%         %some CODE V commands use backslash, which needs to be modified for matlab to write properly
%         backslash = strfind(cmdline,'\');
%         if ~isempty(backslash)
%             cmdline = [cmdline(1:backslash) cmdline(backslash:end)]; %this doubles the backslash \\
%         end
%         
%         fprintf(fid,cmdline);
%         fprintf(fid,'\n');
%     end
%     fclose(fid);


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     