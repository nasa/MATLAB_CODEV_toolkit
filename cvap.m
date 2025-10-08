function [ap_data,lensdata]=cvap
%CVAP gets the aperture data from the current lens
%   
%   function [ap_data,lensdata]=cvap;
%
%   INPUTS: None
%
%   OUTPUT: ap_data: A Nx6 array of raw aperture data
%           lensdata: a strucure lens containing all of the parsed aperture
%           information.
%
%   See also: cvsd

% *****Outputs APERTURE data from CV and imports into MATLAB *************
% This method avoids the CodeV buffer output limitation of requiring the
% output to be numeric. All numeric and string data is passed this way.

MLCVdir = fileparts(mfilename('fullpath'));
cvcmd('buf del b0; BUF YES; APE SA; BUF NO'); % Global Surf Ref = 1
cvcmd('buf sep "|"');
cvcmd(['buf exp sep b0 ' MLCVdir '\temp.txt']);
input = importdata([MLCVdir '\temp.txt'], '\');
delete([MLCVdir '\temp.txt']);
cvcmd('buf fnd b0 "APERTURE DATA"'); %jump to Position line in buffer
istart = cveva('(buf.i)')+2;
cvcmd('buf fnd b0 "BUF NO"'); %jump to IMG line in buffer
iend  = cveva('(buf.i)')-1;

%*****Case for zero aperture data
if length(input{2})>0, disp('No aperture data'); ap_data=0; lensdata=0; return; end 
%*****

% **************** Parses input file into a cell array ********************
for kk=1:iend-5
    rem = char(input(kk+4,:)); 
    for ii=1:3; 
        [token, rem] = strtok(rem,'|');
        [token2,rem2]= strtok(token,' '); 
        if isempty(token2) == 1,
            token2 = ' ';
        end
        if isempty(rem2) == 1,
            rem2 = ' ';
        end
        for jj=1:2
            ap_data(kk,2*ii-1) = cellstr(token2); 
            ap_data(kk,2*ii) = cellstr(rem2); 
end;end;end

for ii=1:size(ap_data,1); 
    if char(ap_data(ii,1)) == char(['ORA']),
        ap_data(ii,2) = ap_data(ii+1,2);
    elseif isempty(char(ap_data(ii,5)))==1,
        ap_data(ii,5) = ap_data(ii,3);
        ap_data(ii,3) = cellstr(['CLR']); % Clear aperture
    elseif isempty(char(ap_data(ii,4)))==1
        ap_data(ii,4) = cellstr([' ']);
    end
    ap_data(ii,2) = strrep(ap_data(ii,2),'S',''); 
end

for jj=1:cvims
    outnum = (find(str2num(char(ap_data(:,2)))==jj));
    for ii=1:size(outnum,1),
        lensdata.surf(jj).ap.raw(ii,1:6) = ap_data(outnum(ii),1:6); 
    end
end
%**************************************************************************

ap_shape = cellstr(['CIR';'REX';'REY';'ELX';'ELY';'ADX';'ADY';'ARO']);
ap_type = cellstr([ 'EDG';'OBS';'HOL']);

% **************** Typce classification for parsing **********************
surfaces=1:cvims;
 for jj=surfaces
     outnum = (find(str2num(char(ap_data(:,2)))==jj));
     if isempty(outnum) == 0,
         for ii=1:size(outnum,1),
             ldat.surf(jj).ap(ii,1:6) = ap_data(outnum(ii),1:6); 
         end
     else  ldat.surf(jj).ap = '0';
     end
end

ap_list = cellstr(['CIR';'REX';'ELX']);
for jj=surfaces
    strout2=[];
    if char(cellstr((ldat.surf(jj).ap(1,1)))) ~='0',
        for ii=1:size(ldat.surf(jj).ap,1);
            adata = char(cellstr(ldat.surf(jj).ap(ii,1)));
            if  (adata==char(ap_list(1)))|...
                (adata==char(ap_list(2)))|...
                (adata==char(ap_list(3))),
                strout2(ii)=1;
            else strout2(ii)=0;
            end
        end
    else strout2=0;
    end
ldat.surf(jj).strout2 = find(strout2);
end

%****************** Aperture separation ***********************************
for jj=surfaces
    for ii=1:size(ldat.surf(jj).strout2,2);
        start = ldat.surf(jj).strout2(ii);
        if ii==size(ldat.surf(jj).strout2,2);
            end1 = size(ldat.surf(jj).ap,1);
        else end1 = ldat.surf(jj).strout2(ii+1)-1;
        end
        ldat.surf(jj).ap(start:end1,6) = num2cell(ii);
    end
end
% *************************************************************************
% Sort aperture data into individual aperertures
for jj=surfaces
for ii=1:size(ldat.surf(jj).strout2,2);
    start = ldat.surf(jj).strout2(ii);
    if ii==size(ldat.surf(jj).strout2,2);
        end1 = size(ldat.surf(jj).ap,1);
    else end1 = ldat.surf(jj).strout2(ii+1)-1;
    end
    if char(ldat.surf(jj).ap(start,1)) == 'REX'
        ldat.surf(jj).mask(ii).type = cellstr('REX');
        ldat.surf(jj).mask(ii).data(:,1) = ldat.surf(jj).ap(start:end1,1);
        ldat.surf(jj).mask(ii).data(:,2) = ldat.surf(jj).ap(start:end1,3);
        ldat.surf(jj).mask(ii).data(:,3) = num2cell(str2num(char(ldat.surf(jj).ap(start:end1,5))));
    elseif char(ldat.surf(jj).ap(start,1)) == 'CIR'
        ldat.surf(jj).mask(ii).type = cellstr('CIR');
        ldat.surf(jj).mask(ii).data(:,1) = ldat.surf(jj).ap(start:end1,1);
        ldat.surf(jj).mask(ii).data(:,2) = ldat.surf(jj).ap(start:end1,3);
        ldat.surf(jj).mask(ii).data(:,3) = num2cell(str2num(char(ldat.surf(jj).ap(start:end1,5))));
    elseif char(ldat.surf(jj).ap(start,1)) == 'ELX'
        ldat.surf(jj).mask(ii).type = cellstr('ELX');
        ldat.surf(jj).mask(ii).data(:,1) = ldat.surf(jj).ap(start:end1,1);
        ldat.surf(jj).mask(ii).data(:,2) = ldat.surf(jj).ap(start:end1,3);
        ldat.surf(jj).mask(ii).data(:,3) = num2cell(str2num(char(ldat.surf(jj).ap(start:end1,5))));
    end
end
end
%**************************************************************************


for ii=surfaces
    lensdata.surf(ii).ap.mask = ldat.surf(ii).mask;
    lensdata.surf(ii).ap.strout = ldat.surf(ii).strout2;
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     