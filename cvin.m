function  [pathfilename,cvoutput] = cvin(pathfilename)
%CVIN  inputs .seq file into CODE V
%
%   function  cvin(pathfilename)
%
%   If pathfilename is left blank, a GUI popup will prompt you to choose a
%   file. If you enter a pathfilename, you do not need to enclose the path
%   in quotes- just enter the path starting with C:\ etc...
%
%   See also CVON, CVOFF, CVCMD, CVSAVE
%

if nargin<1; %no input path to .seq file
    [filename,pathname] = uigetfile({'*.seq', 'SEQ Files (*.seq)';...
                                     '*.len', 'Lens Files (*.len)';...
                                     '*.*', 'All files (*.*)'},'Choose a file'); 
    pathfilename = [pathname filename];
end

end1 = pathfilename(end-2:end);
if (end1(1)=='s'|end1(1)=='S') & (end1(2)=='e'|end1(2)=='E') & (end1(3)=='q'|end1(3)=='Q')
    cmd = ['in "' pathfilename '"'];
     cvoutput = cvcmd(cmd);
elseif (end1(1)=='l'|end1(1)=='L') & (end1(2)=='e'|end1(2)=='E') & (end1(3)=='n'|end1(3)=='N')
    cmd = ['res "' pathfilename '"'];
     cvoutput = cvcmd(cmd);
else
    user_entry = input(sprintf('Do you want to run this file as an SEQ script? \n (y/n)  '),'s');
    if (user_entry(1) == 'y') | (user_entry(1) == 'Y'),
        cmd = ['in "' pathfilename '"'];
        cvoutput = cvcmd(cmd);
    else
        disp(sprintf('Please choose another file. \n'))
    end
end

if nargout<1,
    pathfilename = cvoutput; %allows user to see CodeV output if semicolon is not put after cvin; command
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     