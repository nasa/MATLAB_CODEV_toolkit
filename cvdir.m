function cvdir( dir_option )
% displays directory of CodeV toolbox files
% function cvdir( dir_option )
%   INPUT:  dir_option = 1, CV toolkit m-files
%                      = 2, CV toolkit macro files
%                      = 3, CVUSER folder macro files
%                      = 0, all of the above

if nargin<1, dir_option = 1; end

cvpath = fileparts(which('cvon.m'));

if dir_option == 0 || dir_option == 1
    dir_name = cvpath;
    disp(['Directory: ' dir_name])
    dir(dir_name);
end
if dir_option == 0 || dir_option == 2
    dir_name = [cvpath filesep 'cvmacro'];
    disp(['Directory: ' dir_name])
    dir(dir_name);
end
if dir_option == 0 || dir_option == 3  %FIXME, 
    [~,pathname] = cvsetup;
    if exist(fileparts(pathname),'dir') && exist([fileparts(pathname) filesep 'macros'],'dir')
        dir_name = [fileparts(pathname) filesep 'macros'];
        disp(['Directory: ' dir_name])
        dir(dir_name);
    end
end


