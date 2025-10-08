function pathname = cvpath(functionname)
% CVPATH returns the path of functionname, or of CODE V toolkit if no input
%   function pathname = cvpath(functionname)
%
%   See also CVCMD, CVIN

if nargin<1, functionname = 'cvon'; end

pathname = which(functionname); % Path of functionname called
while pathname(end)~='\', pathname = pathname(1:end-1); end
pathname = pathname(1:end-1); %remove trailing \

pathname2 = fileparts(mfilename('fullpath')); % Path of cvpath

pathname3 = which('cvon'); % Path of cvon
while pathname3(end)~='\', pathname3 = pathname3(1:end-1); end
pathname3 = pathname3(1:end-1); %remove trailing \

if length(pathname) == length(pathname2) & length(pathname) == length(pathname3),
    if pathname == pathname2 & pathname == pathname3,    
    else
        user_entry = input(sprintf(['The returned path does not match the \n'...
            'location of cvpath and cvon \n'...
            'Do you want to continue? \n (y/n)  ']),'s');
        if (user_entry(1) == 'n') | (user_entry(1) == 'N'),
                    pathname=[];
            return
        end
    end
else
    user_entry = input(sprintf(['The returned path does not match the \n'...
            'location of cvpath and cvon \n'...
            'Do you want to continue? \n (y/n)  ']),'s');
    if (user_entry(1) == 'n') | (user_entry(1) == 'N'),
        pathname=[];
        return
    end
end
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     