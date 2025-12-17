function cvdel
%CVDEL deletes CODE V "junk" files that can clutter a directory

delete temp_lens*.len
delete crash_lens*.len
delete *.tmp
delete fort.*
delete cvpma*.int
delete codev*.rec

% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     