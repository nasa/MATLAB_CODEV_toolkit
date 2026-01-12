function [x0,y0,wtf0,vig0,zoom0,xyI0,field0] = cvf(x,y,wtf,vig,zoom,typfld,outstruct)
%CVF retrieves current field points, sets new field points in CodeV, 
%    or sets one (or more) of the current fields already defined
%
%   function [x0,y0,wtf0,vig0,zoom0,xyI0] = cvf(x,y,wtf,vig,z,typfld,outstruct)
%
%   INPUT   x = nan, to return structure of field data
%           x,y = sets new field data, 
%               -- two column input for x is interpreted as x y data
%               -- one column input for x keeps the current field points, 
%                   removing all others not included in the x vector
%           wtf = field weighting
%           vig = vignetting data, size = [numf numz 4]
%           z = zoom position to set data
%           typfld = 'ANG' 'OBJ' 'IMG' or 'RIH'  
%                    or 1 for "ANG", 2 for "OBJ", 3 for "IMG", 4 for "RIH"
%                    [] (default) to use current field type
%                    NOTE: typfld can be input in 2nd or 3rd argument.
%           outstruct = "1" or true to request all field output in a structure
%               ALT INPUT:  cvf(nan) will also provide a structure output.
%
%   OUTPUT  x0 y0 = current field points defined in CodeV (prior to x,y input)
%                  -- if only x ouput is requested, 2 column data x y is given 
%           wtf0 = current field weights
%           vig0 = current vignetting factors for reference rays
%           z0 = zoom position for field data
%           xyI = image coordinate data for field point
%           f0 = field number for data
% 
%   See also CVSD, CVW, CVZ
%


%===========================
if nargin<1, x=[]; end
if nargin<2, y=[]; end
if nargin<3, wtf=[]; end
if nargin<4, vig=[]; end
if nargin<5, zoom = []; end
if nargin<6, typfld = []; end
if nargin<7, outstruct = false; end

N = 14; %numerical precision for CODE V commands

if isstruct(x) %assumes fields structure created by this function
    typfld = x.type;
    wtf = x.Fields_weight;
    vig = x.vig;
    zoom = x.Zoom;
    y = x.y;
    x = x.x; %do x last, since it is the input structure    
end

if isnan(x) 
    outstruct = true; x = []; %simplified input to request structure output
end

if isempty( vig ), vuy=[]; vly=[]; vux=[]; vlx=[]; end
% if isempty( z ), z = 1; end

if ischar(y), typfld = y; y = []; end
if ischar(wtf), typfld = wtf; wtf = []; end

if isempty(typfld), typfld = cveva('(typ fld)',1); end
if ~ischar(typfld)
    if typfld==1, typfld = 'ANG';
    elseif typfld==2, typfld = 'OBJ';
    elseif typfld==3, typfld = 'IMG';
    elseif typfld==4, typfld = 'RIH';
    end
end
        
%retrieve current field data
% tic
% for i_f=1:numf
%     for i_z = 1:numz
%         if contains(typfld,'ANG')
%             x0(i_f) =  cveva(['(xan f' int2str(i_f) 'z' int2str(i_z) ')']);
%             y0(i_f) =  cveva(['(yan f' int2str(i_f) 'z' int2str(i_z) ')']);
%         elseif contains(typfld,'OBJ')
%             x0(i_f) =  cveva(['(xob f' int2str(i_f) 'z' int2str(i_z) ')']);
%             y0(i_f) =  cveva(['(yob f' int2str(i_f) 'z' int2str(i_z) ')']);            
%         elseif contains(typfld,'IMG')
%             x0(i_f) =  cveva(['(xim f' int2str(i_f) 'z' int2str(i_z) ')']);
%             y0(i_f) =  cveva(['(yim f' int2str(i_f) 'z' int2str(i_z) ')']);            
%         elseif contains(typfld,'RIH')
%             x0(i_f) =  cveva(['(xri f' int2str(i_f) 'z' int2str(i_z) ')']);
%             y0(i_f) =  cveva(['(yri f' int2str(i_f) 'z' int2str(i_z) ')']);            
%         end
%     end
% end
% toc

cvmacro('cvfw.seq',1);
[~,numf0,~,numz0] = cvnum;
data = cvbuf(1,1:numf0*numz0,1:11);
x0 = data(:,1);
y0 = data(:,2);
wtf0 = data(:,3);
vuy0 = data(:,4);
vly0 = data(:,5);
vux0 = data(:,6);
vlx0 = data(:,7);
vig0 = [vuy0 vly0 vux0 vlx0];
zoom0 = data(:,8);
xyI0 = data(:,9:10);
field0 = data(:,11);

%format into zooms
x0 =    reshape( x0   , numf0, numz0, []);
y0 =    reshape( y0   , numf0, numz0, []);
wtf0 =  reshape( wtf0 , numf0, numz0, []);
vuy0 =  reshape( vuy0 , numf0, numz0, []);
vly0 =  reshape( vly0 , numf0, numz0, []);
vux0 =  reshape( vux0 , numf0, numz0, []);
vlx0 =  reshape( vlx0 , numf0, numz0, []);
vig0 =  reshape( vig0 , numf0, numz0, []);
zoom0 = reshape( zoom0, numf0, numz0, []);
xyI0 =  reshape( xyI0 , numf0, numz0, []);
field0 =reshape( field0, numf0, numz0, []);

%FIXME:  change format for single zoom output
%                 [numf 4] for single zoom, 4 column data: vuy vly vux vlx

if outstruct
    Fdata.type = cveva('(typ fld)',1);
    Fdata.numF = numf0;
    Fdata.numZ = numz0;
    Fdata.Fields = field0;
    Fdata.Fields_weight = wtf0;
    Fdata.Zoom = zoom0;
    Fdata.x = x0;
    Fdata.y = y0;
    Fdata.fxy = reshape( [x0(:) y0(:)] ,numf0,2,numz0 );
    Fdata.vig = vig0;
    Fdata.Ixy = xyI0;
    Fdata.comments = 'fxy is object space, Ixy is image space';
end

if ~isempty(x) && isempty(y)  %x only entry
    if isscalar(x) || size(x,2)<2 %single column entry meant as field point selection
        f = x;
        x = x0(f);
        y = y0(f);
        wtf = wtf0(f);
        vuy = vuy0(f);
        vly = vly0(f);
        vux = vux0(f);
        vlx = vlx0(f);
        zoom = zoom0(f);
    elseif size(x,2)==2 % parse multi-column x data into y, wtf, vig
%       NOTE:  this needs to be fixed for multiple zooms
        y = x(:,2);
        if size(x,2)>2, wtf = x(:,3); else wtf = []; end
        if size(x,2)>3, vuy = x(:,4); else vuy = []; end
        if size(x,2)>4, vly = x(:,5); else vly = []; end
        if size(x,2)>5, vux = x(:,6); else vux = []; end
        if size(x,2)>6, vlx = x(:,7); else vlx = []; end
        x = x(:,1);
    end
end

if ~isempty(vig)
    vuy = vig(:,:,1);
    vly = vig(:,:,2);
    vux = vig(:,:,3);
    vlx = vig(:,:,4);
end

if ~isempty(x)  %set new field data
    if strcmpi(typfld,'ANG') %object angle
        Xcmd = 'XAN'; Ycmd = 'YAN ';
    elseif strcmpi(typfld,'OBJ') %object height
        Xcmd = 'XOB'; Ycmd = 'YOB ';
    elseif strcmpi(typfld,'IMG') %paraxial image height
        Xcmd = 'XIM'; Ycmd = 'YIM ';
    elseif strcmpi(typfld,'RIH') %real image height
        Xcmd = 'XRI'; Ycmd = 'YRI ';
    else
        disp('ERROR: unknown field type input!  Please use one of: ANG OBJ IMG RIH');
        return
    end
    WTFcmd = []; VUYcmd = []; VLYcmd = []; VUXcmd = []; VLXcmd = [];
    
    numf = length(x);
    numz = 1;
    %set field data (zoom 1)
    for i=1:numf, Xcmd = [Xcmd ' ' num2str(x(i),N)]; end
    cvcmd(Xcmd);
    for i=1:numf, Ycmd = [Ycmd ' ' num2str(y(i),N)]; end
    cvcmd(Ycmd);
    if numz<2 %apply input field data to all zooms, assumes fields aren't zoomed
%         for i=1:length(x), Xcmd = [Xcmd ' ' num2str(x(i),N)]; end
%         cvcmd(Xcmd);
%         for i=1:length(y), Ycmd = [Ycmd ' ' num2str(y(i),N)]; end
%         cvcmd(Ycmd);
        if ~isempty(wtf)
            for i=1:length(wtf), WTFcmd = [WTFcmd 'WTF f' int2str(i) ' ' num2str(wtf(i),N) ';']; end
            cvcmd(WTFcmd);
        end
        if ~isempty(vuy)
            for i=1:length(vuy), VUYcmd = [VUYcmd 'VUY f' int2str(i) ' ' num2str(vuy(i),N) ';']; end
            cvcmd(VUYcmd);
        end
        if ~isempty(vly)
            for i=1:length(vly), VLYcmd = [VLYcmd 'VLY f' int2str(i) ' ' num2str(vly(i),N) ';']; end
            cvcmd(VLYcmd);
        end
        if ~isempty(vux)
            for i=1:length(vux), VUXcmd = [VUXcmd 'VUX f' int2str(i) ' ' num2str(vux(i),N) ';']; end
            cvcmd(VUXcmd);
        end
        if ~isempty(vlx)
            for i=1:length(vlx), VLXcmd = [VLXcmd 'VLX f' int2str(i) ' ' num2str(vlx(i),N) ';']; end
            cvcmd(VLXcmd);
        end
    else  %apply field zoom data  %NEED TO FIX MULTI-ZOOM ENTRY!!!!
        if length(zoom) == 1, zoom = zoom*ones(length(x),1); end
        for i_f=1:numf %number of fields to set
            for i_z=1:numz %number of fields to set
                cmd = [Xcmd ' f' int2str(i_f) ' z' int2str(zoom(i_z)) ' ' num2str(x(i_f),N) ';' Ycmd ' f'  int2str(i_f) ' z' int2str(zoom(i_f)) ' ' num2str(y(i_f),N) ';'];
                cvcmd(cmd);
                if ~isempty(wtf)
                    cmd = [WTFcmd ' f' int2str(i_f) ' z' int2str(zoom(i_z)) ' ' num2str(wtf(i_f),N) ';' Ycmd ' f'  int2str(i_f) ' z' int2str(zoom(i_f)) ' ' num2str(wtf(i_f),N)  ';'];
                    cvcmd(cmd);
                end
                if ~isempty(vuy)
                    cmd = [VUYcmd ' f' int2str(i_f) ' z' int2str(zoom(i_z)) ' ' num2str(vuy(i_f),N) ';' Ycmd ' f'  int2str(i_f) ' z' int2str(zoom(i_f)) ' ' num2str(vuy(i_f),N)  ';'];
                    cvcmd(cmd);
                end
                if ~isempty(vly)
                    cmd = [VLYcmd ' f' int2str(i_f) ' z' int2str(zoom(i_z)) ' ' num2str(vly(i_f),N) ';' Ycmd ' f'  int2str(i_f) ' z' int2str(zoom(i_f)) ' ' num2str(vly(i_f),N)  ';'];
                    cvcmd(cmd);
                end
                if ~isempty(vux)
                    cmd = [VUXcmd ' f' int2str(i_f) ' z' int2str(zoom(i_z)) ' ' num2str(vux(i_f),N) ';' Ycmd ' f'  int2str(i_f) ' z' int2str(zoom(i_f)) ' ' num2str(vux(i_f),N)  ';'];
                    cvcmd(cmd);
                end
                if ~isempty(vlx)
                    cmd = [VLXcmd ' f' int2str(i_f) ' z' int2str(zoom(i_z)) ' ' num2str(vlx(i_f),N) ';' Ycmd ' f'  int2str(i_f) ' z' int2str(zoom(i_f)) ' ' num2str(vlx(i_f),N)  ';'];
                    cvcmd(cmd);
                end
            end
        end
    end
end

% %output only desired zoom data
% x0 =   squeeze( x0   (:,numz,:));
% y0 =   squeeze( y0   (:,numz,:));
% wtf0 = squeeze( wtf0 (:,numz,:));
% vuy0 = squeeze( vuy0 (:,numz,:));
% vly0 = squeeze( vly0 (:,numz,:));
% vux0 = squeeze( vux0 (:,numz,:));
% vlx0 = squeeze( vlx0 (:,numz,:));
% vig0 = squeeze( vig0 (:,numz,:));
% z0 =   squeeze( z0   (:,numz,:));
% xyI0 = squeeze( xyI0 (:,numz,:));
% f0 =   squeeze( f0   (:,numz,:));

if nargout<1 && nargin<1 && length(x0)>1 %plot field data if >1 field point
    for i=1:length(x0(:)), flabel(i) = {[' ' int2str(i)]}; end
    if strcmp(typfld,'ANG')
        figure, scatter(x0(:)*60,y0(:)*60,'.'), text(x0(:)*60,y0(:)*60,flabel), axis equal;
        title('Field defined as Object Angle'); xlabel('Angles in arcminutes');
    else
        figure, scatter(x0(:),y0(:)), text(x0(:),y0(:),flabel), axis equal;
        uni = cvunits;
        if strcmp(typfld,'OBJ'), title('Field defined as Object Height');
        elseif strcmp(typfld,'IMG'), title('Field defined as Paraxial Image Height');
        elseif strcmp(typfld,'RIH'), title('Field defined as Real Image Height');
        end
        if uni==1, xlabel('Field points in mm');
        elseif uni==10, xlabel('Field points in cm');
        else xlabel('Field points in inches');
        end
    end
    x0 = [field0(:) zoom0(:) x0(:) y0(:) wtf0(:) vuy0(:) vly0(:) vux0(:) vlx0(:)];
end

if outstruct
    x0 = Fdata;
elseif nargout==1
    x0 = [x0(:) y0(:)]; 
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     