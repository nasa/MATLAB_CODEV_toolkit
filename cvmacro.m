function [out,cmd] = cvmacro(filename,args)
%CVMACRO  executes an SEQ file with supplied arguments
%
%   function [out,cmd] = cvmacro(filename,args)
%
%   filename = macro name (must be in the \cvmacro directory)
%   args = vector of arguments to be passed to CodeV with the macro
% 
%   See also CVCMD, CVEVA



if nargin<1 
    disp('Must input filename!  cvmacro(filename,args), one of the following:');
    dir([cvpath filesep 'cvmacro']);
    return;
end
if nargin<2, args = {}; end

cmd = ['in "' cvpath '\cvmacro\' filename '" '];
for i=1:length(args)
    if round(args(i))==args(i) %case for integer input
        cmd = [cmd ' ' int2str(args(i))];
    else %use max precision for numerical input
        cmd = [cmd ' ' num2str(args(i),'%20.20f')];
    end
end
out = cvcmd(cmd);

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     