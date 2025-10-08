function [units,unitstring] = cvunits()
%CVUNITS gets units of current lens in CodeV
%
%   function units = cvunits()
%
%   OUTPUT:  units = multiplier from mm, 1.0 10 or 25.4
%            unitstring = 'mm' 'cm' or 'inches'
%        
%   See also CVSD
%

global CodeV

dim = invoke(CodeV,'GetDimension'); % 2=mm, 1=cm, 0=inches
if dim==2, unitstring = 'mm'; units = 1.0;
elseif dim==1, unitstring = 'cm'; units = 10.0;
elseif dim==0, unitstring = 'inches'; units = 25.4;
end


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     