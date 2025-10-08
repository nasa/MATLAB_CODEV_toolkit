function [filename,pathname] = cvsave(pathfilename)
%CVSAVE  saves current lens file under pathfilename
%
%   function [filename,pathname] = CVsave(pathfilename)
%
%   If no input is given, the user will be prompted to select a filename.
%   If pathfilename ends with .seq, the lens will be saved as a script
%   If pathfilename ends with .len, the lens will be saved as a "lens" file
%   If pathfilename ends with neither .seq or .len, two files will be
%       saved, one with the .seq extension and one with the .len extension.
%        
%   See also CVON, CVOFF, CVCMD, CVIN,

if nargin<1, %no input path to .seq file
    [filename,pathname] = uiputfile({'*.seq', 'SEQ Files (*.seq)';...
                                     '*.len', 'Lens Files (*.len)';...
                                     '*.*', 'All files (*.*)'},'Choose a file or write a new filename:');
    pathfilename = [pathname filename];
end

end1 = pathfilename(end-2:end);
if (end1(1)=='s'|end1(1)=='S') & (end1(2)=='e'|end1(2)=='E') & (end1(3)=='q'|end1(3)=='Q')
    cvcmd(['WRL "' pathfilename '"']);
elseif (end1(1)=='l'|end1(1)=='L') & (end1(2)=='e'|end1(2)=='E') & (end1(3)=='n'|end1(3)=='N')
    cvcmd(['SAV "' pathfilename '"']);
else
    savefile = [pathfilename '.seq'];
    cvcmd(['WRL "' savefile '"']);
    savefile = [pathfilename '.len'];
    cvcmd(['SAV "' savefile '"']);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     