function [WL0,WTW0,Wref] = cvw(WL,WTW,outstruct)
%CVW queries and sets wavelengths and weights in CODE V
%
%   function [WL0,WTW0,Wref] = cvw(WL,WTW,outstruct)
%
%   INPUT:  WL, WTW = wavelength data (units in nm) and weights
%            -- two column input for WL is interpreted as WL and WTW
%            -- WL only input sets those current wavelengths, removing
%               all others
%           outstruct = "1" or "true" to output data as a structure
% 
%   OUTPUT: WL0, WTW0 = current wavelengths and weights prior to setting new ones
%           WL0  is a vector of defined wavelengths.
%           WTW0 is a matrix if there are multiple zooms.
%           Wref is a vector of the reference wavelength numbers for each
%                zoom position
%   See also CVSD, CVF, CVZ, CVR

%2011.02.15 NEED TO TEST this for a system with multiple ZOOMs (applying WL
%and WTW

if nargin<1, WL = []; end
if nargin<2, WTW = []; end
if nargin<3, outstruct = false; end

if isnan(WL), outstruct = true; WL = []; end

[~,~,numw,numz] = cvnum;

% ************************** GET WL DATA ********************************
if nargin<2 || nargout>0 || outstruct
    for i_w=1:numw
        WL0(i_w) =  cveva(['(WL w'  int2str(i_w) ')']);
        for i_z = 1:numz
            WTW0(i_w,i_z) = cveva(['(WTW w' int2str(i_w) ' z' int2str(i_z) ')']);
            Wref(i_z) = cveva(['(ref z' int2str(i_z) ')']);
        end
    end
%     tic
%     cvmacro('cvfw.seq',2);
%     data = cvbuf(1,1:numw+1,1:numz+1);
%     WL0 = data(1:numw,1);
%     WTW0 = data(1:numw,2:end);
%     Wref = data(numw+1,2:end);
    if sum(WTW0==0)>0, disp('WARNING: Some wavelength weights are defined as ZERO !'); end
end

if outstruct
    wdata.numW= numw;
    wdata.Wavelengths_nm = WL0;
    wdata.Wavelengths_weight = WTW0;
    wdata.Wref = Wref;
    wdata.zoom_note = 'size weights = [num_W num_Z], size reference = num_Z';
end

% ************************** SET WL DATA ********************************
if ~isempty(WL) % Build CodeV command to set wavelengths
    [~,c] = size(WL);
    if c==2, WTW = WL(:,2); WL = WL(:,1); %convert 2 column input to WL and WTW data
    elseif max(WL)>numw, WL = WL(:); WTW = 0*WL+1; %if number> number of waves, then assume desired wavelength is input
    else WTW = WTW0(WL); WL = WL0(WL); %reset WV values to applicable WV0 values
    end
    if min(WL)<100, disp('Units Check: CodeV expects input wavelengths in units of nm!!!'); end
    WLcmd = 'WL ';
    for i=1:length(WL), WLcmd = [WLcmd ' ' num2str(WL(i))]; end
    cvcmd(WLcmd);

%   CodeV requires integer values for WTW input, so convert if necessary
    changeWTW = 0; 
    for i=1:length(WTW), if ceil(WTW(i))~=WTW(i), changeWTW = 1; end, end
    if changeWTW
        WTW = ceil(WTW*99/max(WTW));
            disp('Converted WTW input to integer data for CodeV');
            disp(WTW);
    end
    WTWcmd = 'WTW ';
    for i=1:length(WTW), WTWcmd = [WTWcmd ' ' num2str(WTW(i))]; end
    cvcmd(WTWcmd);
end

if nargout<1 && nargin<1 && length(WL0)>1 %plot wavelength data
    figure, plot(WL0,WTW0,'o-')
    xlabel('Wavelengths in nm');
    ylabel('Weights');
    title('Spectrum of lens');
    disp(['Wavelengths and Weights']);
    disp([WL0(:) WTW0(:)]);
    WL0 = [];
end

if outstruct
    WL0 = wdata;
end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     