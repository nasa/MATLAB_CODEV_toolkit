function [data,inputs] = cvr(arg1,arg2,arg3,arg4,arg5,arg6)
%CVR  gets CodeV ray trace data for reference rays from database
%
% FOR REFERENCE RAY DATA
%   function data = cvr(raydata,f,r,s,w,z)
%   
%   INPUT V1:  raydata = text input X,Y,Z,L,M,N,AOI,AOR,SRL,SRM,SRN,OPL,or OPD 
%              field, default = 1
%              ray, 1,2,3,4, or 5:  chief,top,bottom,left,right, default = 1
%              surface, default = sI (image), can be multiple surfaces (vector input)
%              wave, default = 1
%              zoom, default = 1
%   OUTPUT V1: data for each surface input
%
%   INPUT V2: vector of [f r s w z], default = [1 1 image 1 1]
%   OUTPUT V2: vector of [X,Y,Z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP]
%               OP from s1..s
%
% FOR RELATIVE PUPIL OR FIELD
%   function data = cvr(f,r,s,w,z)
%
%   INPUT:  f = Nx2 vector of relative field coordinates. If you use [] to leave
%               the field specification blank, the default is [0,0].
%           r = Nx2 vector of realtive pupil coordinates. If you use [] to leave
%               the pupil specification blank, the default is [0,0].
%           s = vector of surfaces to be traced.
%           w = vector of wavelegnths to be traced.
%           z = vector of zooms to be traced.
%       NOTE: Only one input of s, w, or z, can be a vector per run. If f/r
%       are larger than 1x2 and s,w,or z is a vector, the output will be a
%       cell array with each element containing all f/r data for a given
%       s,w,z parameter.
%
%   OUTPUT: 
%   X,Y,Z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP
%
%
%   See also CVRGRID
%


datatype = {'x','y','z','l','m','n','aoi','aor','srl','srm','srn','op'};
[gets,getf,getw,getz] = cvnum;

if nargin<1, arg1 = 0; end 

    if nargin<1, f=1; r=1; s=cvnum; w=1; z=1; %*************** NO ARGS ********
    elseif nargin==1 | (nargin==2 & ischar(arg2)), %******** ONE ARG ********* 
        raydata = arg1;
        if size(raydata,2)>0, f = raydata(:,1); else f = 1; end
        if size(raydata,2)>1, r = raydata(:,2); else r = 1; end
        if size(raydata,2)>2, s = raydata(:,3); else s = cvnum; end
        if size(raydata,2)>3, w = raydata(:,4); else w = 1; end
        if size(raydata,2)>4, z = raydata(:,5); else z = 1; end 
    elseif nargin==5 | (nargin==6 & ischar(arg6)) %******** FIVE ARGS ********
        if isempty(arg1), f=1:getf; else f = arg1; end
        if isempty(arg2), r=1:5; else r = arg2; end
        if isempty(arg3), s=1:gets; else s = arg3; end
        if isempty(arg4), w=1:getw; else w = arg4; end
        if isempty(arg5), z=1:getz; else z = arg5; end
    end
    
    [rownum,c]=size(f);
    if rownum>1 & c>1, %relative field coordinates
        relf = 'y';
        if rownum~=2 & c==2, 
            f=f'; 
        elseif rownum==2 & c==2, 
            disp('NOTE: make sure relative field coordinates are along the second dimension'); 
        end
        numf = rownum;
    else   
        relf='n';
        numf=length(f); 
    end
    
    [rownum,c]=size(r);
    if rownum>1 & c>1, %relative pupil coordinates
        relpup = 'y';
        if rownum~=2 & c==2, 
            r=r';
        elseif rownum==2 & c==2, 
            disp('NOTE: make sure relative pupil coordinates are along the second dimension'); 
        end
        numr = rownum;
    else
        relpup='n';
        numr=length(r);
    end
    
    numr=length(r); nums=length(s); numw=length(w); numz=length(z);
    max1 = max([numf numr nums numw numz]);
    f = padarray_blu(f,[max1-numf]);
    r = padarray_blu(r,[max1-numr]);
    s = padarray_blu(s,[max1-nums]);
    w = padarray_blu(w,[max1-numw]);
    z = padarray_blu(z,[max1-numz]);
    input = [f;r;s;w;z]';
    
    cvbufput(input,8);
    rows1 = numf*numr*nums*numw*numz;
    if relf~='y',
        if relpup~='y',
            cvmacro('cvr.seq', [1 numf numr nums numw numz 8 1]); %first argument specifies type
        elseif relpup=='y', %relative pupil coordinates
            cvmacro('cvr.seq', [3 numf numr nums numw numz 8 1]); %first argument specifies type
        end
    elseif relf=='y', %relative field coordinates
        if relpup~='y',
            cvmacro('cvr.seq', [2 numf numr nums numw numz 8 1]); %first argument specifies type
        elseif relpup=='y', %relative pupil coordinates
            cvmacro('cvr.seq', [4 numf numr nums numw numz 8 1]); %first argument specifies type
        end
    end
    output = cvbuf(1,1:rows1,1:12);

    ii=1;
    for ff=1:numf; for rr=1:numr; for ss=1:nums; for ww=1:numw; for zz=1:numz
                        data(:,ff,rr,ss,ww,zz) = output(ii,:); %Max is a 6D array containing all data
                        ii=ii+1;
                    end; end; end; end; end

    %**************** If only one datatype is specified *********
    if nargin==2 & ischar(arg2),
        for ii=1:12
            if  arg2(1)==datatype{ii}(1) & arg2==datatype{ii}(1:length(arg2)),
                data = data(ii,:,:,:,:,:);
            end
        end
    elseif nargin==6 & ischar(arg6),
        for ii=1:12
            if arg6(1)==datatype{ii}(1) & arg6==datatype{ii}(1:length(arg6)),
                data = data(ii,:,:,:,:,:);
            end
        end
    end
    data = squeeze(squeeze(squeeze(data))); % squeeze out all singlton dimensions

    inputs.f=f;
    inputs.r=r;
    inputs.s=s;
    inputs.w=w;
    inputs.z=z;
end
        
function out = padarray_blu(arg1,arg2);
    out = ones(1,length(arg1)+arg2);
    out(1:length(arg1)) = arg1;
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     