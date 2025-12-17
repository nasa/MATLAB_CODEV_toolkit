function [filename,pathname,out] = cvsave(pathfilename,seq_strings)
%CVSAVE  saves current lens file under pathfilename
%
%   function [filename,pathname] = CVsave(pathfilename,seq_strings)
%   INPUT:  pathfilename, 
%               If no input is given, the user will be prompted to select a filename.
%               If pathfilename ends with .seq, the lens will be saved as a script
%               If pathfilename ends with .len, the lens will be saved as a "lens" file
%               If pathfilename ends with neither .seq or .len, two files will be
%                   saved, one with the .seq extension and one with the .len extension.             
%           seq_strings = string to save as a SEQ file, 
%                       OR
%                 {file1,file2,...} to concatenate into a single SEQ file
%   See also CVON, CVOFF, CVCMD, CVIN,

if nargin<1 %no input path to .seq file
    [filename,pathname] = uiputfile({'*.seq', 'SEQ Files (*.seq)';...
                                     '*.len', 'Lens Files (*.len)';...
                                     '*.*', 'All files (*.*)'},'Choose a file or write a new filename:');
    pathfilename = [pathname filename];
end
if nargin<2, seq_strings = []; end


if isempty(seq_strings)
    [~,~,EXT] = fileparts(pathfilename);
    if contains(upper(EXT),'SEQ')
        out = cvcmd(['WRL "' pathfilename '"']);
    elseif contains(upper(EXT),'LEN')
        out = cvcmd(['SAV "' pathfilename '"']);
    else
        savefile = [pathfilename '.len'];
        cvcmd(['SAV "' savefile '"']);
        savefile = [pathfilename '.seq'];
        out = cvcmd(['WRL "' savefile '"']);
    end
else
    if iscell(seq_strings)
        seq_files = seq_strings;
        seq_strings = [];
        for i=1:length(seq_files)
            seq_strings = [seq_strings fileread(seq_files{i})];
        end
    end
    writelines(seq_strings,pathfilename);
end


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     