function cvexp2ims(data)
% CVEXP2IMS converts pupil surface to image
% function data = cvexp2ims(data)

if nargin<1, disp('Please include data: cvexp2ims(data)'); end

cvcmd(['RDY SI ' num2str(data(10),20)]);
cvdec(cvnum, data(1:9));

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     