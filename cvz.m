function activezooms = cvz(z)
%CVZ queries and sets active zoom positions in CODE V
%
%   function activezooms = cvz(z)
%
%   INPUT:  z = vector of desired active zoom positions
%   OUTPUT: activezooms = current active zoom positions
%           
%   See also CVSD, CVF, CVW
%

[nums,numf,numw,numz] = cvnum;

cvcmd('buf no; buf del b0; BUF YES; POS ?; BUF NO'); % Output to buffer b0
onoff = cveva('(buf.txt b0 i3)',1); %get buffer line with zoom data
activezooms = [];
for i=1:numz
    [test, onoff] = strtok(onoff);
    if strcmp(test,'ON'), activezooms = [activezooms i]; end
end

%CODEV syntax note
% to activate zoom 1 and turn zoom 2 off:  pos z1 on; pos z2 off; 
% one zoom must be on at all times, or codev will reject the command

if nargin>0
    Zcmd = ['POS ']; %currently sets all zoom positions active
    if length(find(z>numz))>0 | length(find(z<1)>0), disp('Requested Zoom Position does not exist.'), return, end;
    z = sort(z); %ensure values are in ascending order
    for j=1:numz
        if length(find(z==j))>0, Zcmd = [Zcmd ' y'];
        else Zcmd = [Zcmd ' n'];
        end
    end
    cvcmd(Zcmd);
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     