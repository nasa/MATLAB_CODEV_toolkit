%CVOFF kills the COM link between CODE V and Matlab
%
%
%   See also CVON, CVCMD, CVDB, CVIN, CVOPEN, CVSAVE 
%

%note: this was a function, but is now a script

if exist('CodeV','var') 
    disp('Stopping CodeV...');
%     invoke(CodeV,'StopCodeV'); %pre v11.3
    CodeV.release; %I think this removes the COM link
    clearvars -GLOBAL CodeV;
    clear CodeV;
    disp('CodeV stopped.');
else
    clearvars -GLOBAL CodeV;
    clear CodeV;
end
% evalin('base','clear global CodeV'); %pre v11.3

if exist('CODEVserver','var')
    clear CODEVserver;
end

cvdel; %clear junk files

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22      