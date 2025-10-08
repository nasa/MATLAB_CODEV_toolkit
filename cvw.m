function [WL0,WTW0] = cvw(WL,WTW)
%CVW queries and sets wavelengths and weights in CODE V
%
%   function [WL0,WTW0] = cvw(WL,WTW)
%
%   INPUT:  WL, WTW = wavelength data (units in nm) and weights
%            -- two column input for WL is interpreted as WL and WTW
%            -- WL only input sets those current wavelengths, removing
%               all others
%   OUTPUT: WL0, WTW0 = current wavelengths and weights prior to setting new ones
%
%   See also CVSD, CVF, CVZ, CVR

[nums,numf,numw] = cvnum;

% ************************** GET WL DATA ********************************
if nargin<2 || nargout>0
    cvmacro('cvfw.seq',2);
    data = cvbuf(1,1:numw,1:2);
    WL0 = data(:,1);
    WTW0 = data(:,2);
end

% ************************** SET WL DATA ********************************
if nargin>0 % Build CodeV command to set wavelengths
    if nargin==1, 
        [r,c] = size(WL);
        if c==2, WTW = WL(:,2); WL = WL(:,1); %convert 2 column input to WL and WTW data
        elseif max(WL)>numw, WL = WL(:); WTW = 0*WL+1; %if number> number of waves, then assume desired wavelength is input
        else WTW = WTW0(WL); WL = WL0(WL); %reset WV values to applicable WV0 values
        end 
    end
    if min(WL)<100, disp('Units Check: CodeV expects input wavelengths in units of nm!!!'); end
    WLcmd = ['WL '];
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
    WTWcmd = ['WTW '];
    for i=1:length(WTW), WTWcmd = [WTWcmd ' ' num2str(WTW(i))]; end
    cvcmd(WTWcmd);
end

if nargout<1 && nargin<1 %plot wavelength data
    figure, scatter(WL0,WTW0)
    xlabel('Wavelengths in nm');
    ylabel('Weights');
    title('Spectrum of lens');
    disp(['Wavelengths and Weights']);
    disp([WL0' WTW0']);
end

if nargout==1, WL0 = [WL0(:) WTW0(:)]; end

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     