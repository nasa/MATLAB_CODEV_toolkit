function [x0,y0] = cvf(x,y)
%CVF retrieves current field points, sets new field points in CodeV, 
%    or sets one (or more) of the current fields already defined
%
%   function [x0,y0] = cvf(x,y)
%
%   INPUT x,y = sets new field data, 
%               -- two column input for x is interpreted as x y data
%               -- x only input sets those current field points, removing
%               all others
%
%   OUTPUT x0 y0 = current field points defined in CodeV (prior to x,y input)
%                  -- if only x ouput is requested, 2 column data x y is given 
%
%   See also CVSD, CVW, CVZ
%

%===========================
typfld = cveva('(typ fld)',1);
[nums,numf] = cvnum;

if nargin<2 || nargout>0 %retrieve current field data if needed
    cvmacro('cvfw.seq',1);
    data = cvbuf(1,1:numf,1:2);
    x0 = data(:,1);
    y0 = data(:,2);
end

if nargin>0      % Build CodeV command to set fields
    if nargin==1, 
        [r,c] = size(x);
        if c==2, y = x(:,2); x = x(:,1); %convert 2 column input to x y data
        else y=y0(x); x=x0(x); %set existing field points
        end
    end

    if strcmp(typfld,'ANG') %object angle
        Xcmd = ['XAN']; Ycmd = ['YAN '];
    elseif strcmp(typfld,'OBJ') %object height
        Xcmd = ['XOB']; Ycmd = ['YOB '];
    elseif strcmp(typfld,'IMG') %paraxial image height
        Xcmd = ['XIM']; Ycmd = ['YIM '];
    elseif strcmp(typfld,'RIH') %real image height
        Xcmd = ['XRI']; Ycmd = ['YRI '];
    end
    
    for i=1:length(x), Xcmd = [Xcmd ' ' num2str(x(i),20)]; end
    cvcmd(Xcmd);
    for i=1:length(y), Ycmd = [Ycmd ' ' num2str(y(i),20)]; end
    cvcmd(Ycmd);
    
end

if nargout<1 & nargin<1 %plot field data
    for i=1:length(x0), flabel(i) = {[' ' int2str(i)]}; end
    if strcmp(typfld,'ANG')
        figure, scatter(x0*60,y0*60,'.'), text(x0*60,y0*60,flabel), axis equal;
        title('Field defined as Object Angle'); xlabel('Angles in arcminutes');
    else
        figure, scatter(x0,y0), text(x0,y0,flabel), axis equal;
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
end

if nargout==1, x0 = [x0' y0']; end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     