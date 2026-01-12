function [activezooms,titles] = cvz(z,outstruct)
%CVZ queries and sets active zoom positions in CODE V
%
%   function [activezooms,titles] = cvz(z,outstruct)
%
%   INPUT:  z = vector of desired active zoom positions
%           outstruct = true to output structure format
%   OUTPUT: activezooms = current active zoom positions
%           titles = names of zoom positions
% 
%   See also CVSD, CVF, CVW
%

[~,~,~,numz] = cvnum;
if nargin<2, outstruct = false; end

% if nargout<1 || nargout>1
    for i_z = 1:numz
        z_str = int2str(i_z);
        activezooms(i_z) = cvdb(['pos z' z_str]);
        titles{i_z} = cvdb(['tit z' z_str],1);
    end
% end

if outstruct
    Zdata.numZ = numz;
    Zdata.Zoom_Active = activezooms;
    Zdata.Zoom_Titles = titles;
end

%CODEV syntax note
% to activate zoom 1 and turn zoom 2 off:  pos z1 on; pos z2 off; 
% one zoom must be on at all times, or codev will reject the command

if nargin>0
    Zcmd = ['POS ']; %currently sets all zoom positions active
    if ~isempty(find(z>numz, 1)) || ~isempty(find(z<1)>0), disp('Requested Zoom Position does not exist.'), return, end;
    z = sort(z); %ensure values are in ascending order
    for j=1:numz
        if ~isempty(find(z==j, 1)), Zcmd = [Zcmd ' y'];
        else Zcmd = [Zcmd ' n'];
        end
    end
    cvcmd(Zcmd);
end

if nargout<1 && nargin < 1
    disp(' ');
    disp('z  on   title');
    disp('=============');
    for i_z = 1:numz
        txt = [int2str(i_z) '   ' int2str(activezooms(i_z)) '   ' titles{i_z}];
        disp(txt);
    end
end

if outstruct, activezooms = Zdata; end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     

% ARCHIVE CODE

% cvmacro('cvfw.seq',3);
% if nargout<1 || nargout>1
%     data_titles = cvbuf(1); %gets entire buffer
%     if numz<2
%         data = textscan(data_titles,'%d %s','Delimiter','');
%         znum = data{1};
%         titles = data{2};
%     else
%         data = cvbuf(1,1:numz,1:2);
%         titles = cvbufstr(1,1:numz,3);
%         znum = data(:,1);
%         on_off = data(:,2);
% %         for i_z = 1:numz
% %             [znum(i_z),on_off(i_z),titles{i_z}] = strread(data_titles(i_z,:),'%d %d %s');
% %         end
%     end
% end
% if nargout>0
%     activezooms = cvbuf(1,1:numz,1);
% end


% cvcmd('buf no; buf del b0; BUF YES; POS ?; BUF NO'); % Output to buffer b0
% onoff = cveva('(buf.txt b0 i3)',1); %get buffer line with zoom data
% activezooms = [];
% for i=1:numz
%     [test, onoff] = strtok(onoff);
%     if strcmp(test,'ON'), activezooms = [activezooms i]; end
% end

% if nargout>1 || nargout<1
%     for i=1:numz
%         titles{i} = cveva(['(tit z' int2str(i) ')'],1);
%     end
% end
% if nargout<1
%     for i=1:numz
%         disp([int2str(i) '   ' titles{i}]);
%     end
% end
