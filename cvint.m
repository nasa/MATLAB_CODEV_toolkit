function [output1,output2] = cvint(arg1,arg2,arg3,arg4,arg5)
%CVINT  attaches an INT file to optical surface or reads in an INT file to Matlab
%
%option 1:   function [data,mask] = cvint(s,pathfilename,options)
%option 2:   function pathfilename = cvint(s,data,mask,options,pathfilename)
%option 3:   function [data,mask] = cvint(pathfilename)
%
%   INPUTS: s = surface on which to apply INT file
%           pathfilename = name of INT file to be read/written (if left
%               blank, the function will prompt you to choose a file
%           options = a cell array (i.e. use curly brackets {} ) containing
%               choices for applying INT to surface
%               options = {type,label,radius,scale,decenter_x,decenter_y,rotation,mirrordata};
%               (default: {'sur',filename,cvmap,1,0,0,0,'NO'} )
%                   type: SUR, WRF, FIL
%           array = if argument #2 is numeric, the array will be converted
%               to an INT file and then applied to a surface.
%               *** array is in units of meters and is applied as a surface map 
%           mask = mask for INT file
%
%           NOTE: If options is left empty (i.e. use {} ), options will be
%           defaults. If you use [] for pathfilename in option 2, INT file
%           names will automatically be generated. Also, option 2 works for
%           vector surface inputs and 3D array inputs.
%
%   OUTPUT: pathfilename = path and file name of INT file read or saved
%           data = when reading an INT file, reaturns numeric data
%           mask = array of zeros and ones corresponding to actual data
%
%   See also:

if nargin<1 | ischar(arg1), % Read an INT file into Matlab   
    if nargin<1,
        [filename,pathname] = uigetfile({'*.int',  'int Files (*.int)'},'Pick an .int file');
        pathfilename = [pathname filename];
    else 
        pathfilename = arg1;
        [pathname,a,b]=fileparts(arg1);
        filename = [a b];
    end
    fileidread = fopen(pathfilename,'rt');
    inputstr = fgets(fileidread);
    while inputstr(1)=='!', inputstr = fgets(fileidread); end
    INTtitle = inputstr;
    INTformat = fgets(fileidread); %Format line for .int file
    intdata = fscanf(fileidread,'%d',inf); % .int data itself
    fclose(fileidread);
    cellDATA = strread(INTformat,'%s','delimiter',' '); %chop up input into cell
    for i=1:length(cellDATA),
        if strfind(cellDATA{i},'GRD')==1,  gridsize = str2num(cellDATA{i+1}); gridsizeY = str2num(cellDATA{i+2}); end
        if strfind(cellDATA{i},'WVL')==1,  wavelength = str2num(cellDATA{i+1}); end
        if strfind(cellDATA{i},'NDA')==1,  NoDataValue = str2num(cellDATA{i+1}); end
        if strfind(cellDATA{i},'SSZ')==1,  scalesize = str2num(cellDATA{i+1}); end
        if strfind(cellDATA{i},'XSC')==1,  xscale = str2num(cellDATA{i+1}); end
    end
    intdata = reshape(intdata,gridsize,gridsizeY);
%       transpose data for proper orientation
    intdata(:,:) = flipud(intdata(:,:)');
    data = find(intdata~=NoDataValue);
    nodata = find(intdata==NoDataValue);
    intdata(nodata) = 0; %set NoDataValue to zero
    mask=ones(size(intdata,1),size(intdata,2));
    mask(nodata)=0;

%       convert int data to units of meters
    intdata = intdata*(wavelength/scalesize)/1e6;
    if nargout<1
        assignin('caller',filename(1:end-4),intdata);
        assignin('caller',[filename(1:end-4) '_mask'],mask);
    end
    output1 = intdata;
    output2 = mask;
    
%**************************************************************************
%******************** Apply INTs to surface *******************************
else
    if size(arg2,3)>length(arg1),
        s=ones(size(arg2,3),1); arg1=arg1.*s;
    elseif size(arg2,3)<length(arg1),
        for kk=1:length(arg1),
            temp(:,:,kk)=arg2(:,:,1);
        end
        arg2=temp;
    end
for ii=1:length(arg1),
    
    s = arg1(ii);
    if nargin<2 | isempty(arg2), % Applies an INT file to a surface
        [filename,pathname] = uigetfile({'*.int',  'int Files (*.int)'},'Choose an .int file'); 
        pathfilename = [pathname filename];
        output1 = pathfilename;
        if nargin<3, % Options for INT read
            options = {'sur',filename,cvmap,1,0,0,0,'NO'};
        elseif nargin>3 & ~iscell(arg4),
            disp('"options" must be a cell array');
            return
        elseif iscell(arg3) | iscell(arg4),
            if nargin==3, options=arg3; elseif nargin>3, options=arg4; end
            if length(options)<1, options{1} = 'sur'; end
            if length(options)<2, options{2} = filename; end
            if length(options)<3, options{3} = cvmap; end
            if length(options)<4, options{4} = 1; end
            if length(options)<5, options{5} = 0; end
            if length(options)<6, options{6} = 0; end
            if length(options)<7, options{7} = 0; end
            if length(options)<8, options{8} = 'NO'; end
        end
        if nargout>0, [output1,output2]=cvint(pathfilename);end
    elseif ischar(arg2), % Applies an INT file to a surface
        pathfilename = arg2;
        [pathname,a,b]=fileparts(arg2);
        filename = [a b];
        output1 = pathfilename;
        if nargin<3, % Options for INT read
            options = {'sur',filename,cvmap,1,0,0,0,'NO'};
        elseif nargin>3 & ~iscell(arg4),
            disp('"options" must be a cell array');
            return
        elseif iscell(arg3) | iscell(arg4),
            if nargin==3, options=arg3; elseif nargin>3, options=arg4; end
            if length(options)<1, options{1} = 'sur'; end
            if length(options)<2, options{2} = filename; end
            if length(options)<3, options{3} = cvmap; end
            if length(options)<4, options{4} = 1; end
            if length(options)<5, options{5} = 0; end
            if length(options)<6, options{6} = 0; end
            if length(options)<7, options{7} = 0; end
            if length(options)<8, options{8} = 'NO'; end
        end
        if nargout>0, [output1,output2]=cvint(pathfilename);end
%**************************************************************************
%**************************************************************************
    elseif isnumeric(arg2), % creates an INT file from an array and applies it to a surface
        if nargin<5,  % Creates a new file
            [filename,pathname] = uiputfile({'*.int',  'int Files (*.int)'},'Choose a name for the .int file');
            end1 = filename(end-2:end);
            if ~((end1(1)=='i'|end1(1)=='I') & (end1(2)=='n'|end1(2)=='N') & (end1(3)=='t'|end1(3)=='T')),
                filename = [filename '.int'];
            end
            pathfilename = [pathname filename];
        elseif ischar(arg5), % Applies an INT file to a surface
            pathfilename = arg5;
            [pathname,a,b]=fileparts(arg5);
            filename = [a b];
        elseif isempty(arg5),
            pathname = 'c:\cvuser\';
            filename = ['S' num2str(s) '_INT_' num2str(ii) '.int'];
            pathfilename = [pathname filename];
        else
            disp('Did not complete successfully');
        end
        if nargin<4, % Options for INT read
            options = {'sur',filename,cvmap,1,0,0,0,'NO'};
        elseif nargin>3 & ~iscell(arg4),
            disp('"options" must be a cell array');
            return
        elseif iscell(arg3) | iscell(arg4),
            if nargin==3, options=arg3; elseif nargin>3, options=arg4; end
            if length(options)<1, options{1} = 'sur'; end
            if length(options)<2, options{2} = filename; end
            if length(options)<3, options{3} = cvmap; end
            if length(options)<4, options{4} = 1; end
            if length(options)<5, options{5} = 0; end
            if length(options)<6, options{6} = 0; end
            if length(options)<7, options{7} = 0; end
            if length(options)<8, options{8} = 'NO'; end
        end
        
        [w,wt] = cvw;
        w = w.*1e-9; % wavelength in meters
        inttype = options{1};
        title = pathfilename;
        comments = 'int file generated using parts of makeintfile.m function';
        data = arg2(:,:,ii);
        mask = arg3;
%        mask = ones(size(data,1),size(data,2)); mask(data==0)=0;

% NOTE: the input array data is assumed to be generated by matlab so that it
% is read from -x to +x along a row, and -y to + y from top to bottom
% of a column.  So it needs to be transposed and flipped before mapped to 1D array
        data = fliplr(transpose(data));%for correctly mapped in 1D array
        mask = fliplr(transpose(mask));
        nodataval = 32767; %default CodeV value for nodata
        [nx,ny]=size(data);
        maxdata = max(max(abs(data(mask==1))));
        two15m2 = 2^15 - 2; %32766
        if maxdata < two15m2*w % maximum spread of data where scalesize = 1
            data16 = floor(two15m2*data/maxdata).*mask + nodataval.*(~mask); % minimum step size to fill range of .int data
            ssz = floor(two15m2*w/maxdata);
        else
            disp('Data exceeds maximum amplitude of (2^15 - 2) * lambda !'), return,
        end

        fid = fopen(pathfilename,'w'); fprintf(fid,'!Precision of this .int file = %f m\n',ssz*maxdata/two15m2);
        fprintf(fid,'!%s \n',comments); fprintf(fid,'%s \n',title);
        fprintf(fid,['GRD %d %d %s WVL %10.5f SSZ %6d  NDA ' int2str(nodataval) ' \n'],nx,ny,inttype,w*1e6,ssz);
        fprintf(fid,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d \n',data16);
        fclose(fid);
        if length(arg1)>1,
            output1{ii}=pathfilename;
        else
            output1=pathfilename;
        end
    end



    
%**************************************************************************        

    cvcmd([...
        'INT S' num2str(s) ' L''' options{1} ''' ' pathfilename ';' ]);
    cvcmd([...
        'ISF S' num2str(s) ' L''' options{1} ''' ' num2str(options{4}) ';' ...
        'INR S' num2str(s) ' L''' options{1} ''' ' num2str(options{3}) ';' ...
        'INX S' num2str(s) ' L''' options{1} ''' ' num2str(options{5}) ';' ...
        'INY S' num2str(s) ' L''' options{1} ''' ' num2str(options{6}) ';' ...
        'IRO S' num2str(s) ' L''' options{1} ''' ' num2str(options{7}) ';'...
        'IMI S' num2str(s) ' L''' options{1} ''' ' options{8} ';' ...
        ]);
end
end


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     