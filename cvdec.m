function [dec,typ] = cvdec(s,setval)
%CVDEC  returns decenter data from CodeV
%
%   function [dec,typ] = cvdec(s,setval);
%
%   INPUTS: s = surfaces desired, default = all, including object
%           setvals =  matrix of values to be shifted. You only need
%            to specify out as many columns are there are data. I.E.: For
%            only x- and y- shifts, you only need 2 columns: [XDE; YDE]
%            Columns are in this order: 
%            1=XDE, 2=YDE, 3=ZDE, 4=ADE, 5=BDE, 6=CDE, 7=XOD, 8=YOD, 9=ZOD
%
%   OUTPUT: dec = array, row data = surfaces
%                 column data = xde yde zde ade bde cde xod yod zod
%   NOTE: Object data is returned in first row when no input s vector is provided.
%         Thus, surface 1 data will be row 2, s3 = row 4, etc.
%           
%   See also, CVSENS, CVLOM, CVRAYGRID

nums = cvnum;
if nargin<1, s=1:nums; end

cvmacro('cvdec.seq');
dec = cvbuf(1,1:nums+1,1:9);

if nargout~=1 %get typ of surface decenter
   for i=1:nums+1
       typ{i} = cveva(['(buf.str b1 i' int2str(i) ' j10)'],1);
   end
   typ = typ(s+1); %return only surfaces desired
end
if nargin>0
    dec = dec(s+1,:); %return only surfaces desired
end
%**************************************************************************
if nargin==2,
    [r,c] = size(setval);
    XDE = setval(:,1);
    if c<2, YDE = dec(:,2); else YDE = setval(:,2); end;
    if c<3, ZDE = dec(:,3); else ZDE = setval(:,3); end;
    if c<4, ADE = dec(:,4); else ADE = setval(:,4); end;
    if c<5, BDE = dec(:,5); else BDE = setval(:,5); end;
    if c<6, CDE = dec(:,6); else CDE = setval(:,6); end;
    if c<7, XOD = dec(:,7); else XOD = setval(:,7); end;
    if c<8, YOD = dec(:,8); else YOD = setval(:,8); end;
    if c<9, ZOD = dec(:,9); else ZOD = setval(:,9); end;

    for ii=1:r
        sn = num2str(s(ii));
        if (XOD(ii)==0) & (YOD(ii)==0) & (ZOD(ii)==0) ,
            ODstr = ' ';
        else ODstr = ['xod s' sn ' ' num2str(XOD(ii),20) ...
                '; yod s' sn ' ' num2str(YOD(ii),20) ...
                '; zod s' sn ' ' num2str(ZOD(ii),20) ];
        end
        cmd = ['xde s' sn ' ' num2str(XDE(ii),20) ...
            '; yde s' sn ' ' num2str(YDE(ii),20) ...
            '; zde s' sn ' ' num2str(ZDE(ii),20) ...
            '; ade s' sn ' ' num2str(ADE(ii),20) ...
            '; bde s' sn ' ' num2str(BDE(ii),20) ...
            '; cde s' sn ' ' num2str(CDE(ii),20) ...
            ';' ODstr ];
        cvcmd(cmd);
    end

end

if nargout<1
    for i=1:length(typ)
        data{i,1} = typ{i};
        for j=1:9
            data{i,j+1} = dec(i,j);
        end
    end
    dec = data
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     