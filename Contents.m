% CODEV
%
% Files
%   cvap       - gets the aperture data from the current lens
%   cvbestsph  - returns radius of best fit sphere
%   cvbpr      - gets CodeV beam propagation data
%   cvbuf      - gets buffer contents from CODE V
%   cvbufput   - puts array data into the CodeV buffer
%   cvcd       - changes current directory in matlab to CodeV toolkit directory
%   cvcmd      - sends command to CODE V command line
%   cvcv       - starts CODE V GUI version and loads file
%   cvdb       - queries CODE V database item
%   cvdec      - returns decenter data from CodeV
%   cvdraw     - draws the current lens in Matlab
%   cvenc      - gets CODE V PSF based encircled energy for fieldpoint 1
%   cveva      - sends command to CODE V command line over existing COM link
%   cvexp      - uses CodeV RAYTRA engine to calculate ray data at exit pupil
%   cvexp2ims  - converts pupil surface to image
%   cvf        - retrieves current field points, sets new field points in CodeV, 
%   cvfl       - evaluates the focal length
%   cvgc       - gets the global coodinates for surfaces, and converts using makeglbs.seq
%   cvhelp     - displays tookit functions by subject
%   cvims2exp  - converts image surface to exit pupil
%   cvin       - inputs .seq file into CODE V
%   cvint      - attaches an INT file to optical surface or reads in an INT file to Matlab
%   cvlicense  - shows the NASA OPEN SOURCE SOFTWARE AGREEMENT
%   cvmacro    - executes an SEQ file with supplied arguments
%   cvmap      - returns max aperture of surface s, zoom z, in CodeV
%   cvnum      - gets number of surfaces, field points, wavelengths, and zooms 
%   cvoff      - kills the COM link between CODE V and Matlab
%   cvon       - starts the COM link between Matlab and CODE V
%   cvopl      - uses CodeV RAYTRA engine to calculate OPL data at image surface.
%   cvout      - retrieves CodeV output
%   cvpath     - returns the path of functionname, or of CODE V toolkit if no input
%   cvpin      - gets CodeV exit pupil intensity data
%   cvpma      - gets CODE V exit pupil map, mask, and other PMA data
%   cvpsf      - gets CodeV PSF normalized to perfect lens (i.e. strehl value)
%   cvr        - gets CodeV ray trace data for reference rays from database
%   cvrac      - gets ray failure map from CODE V
%   cvrbshift  - "shifts" the decenters and tilts for single surface 'surfnum', 
%   cvreadme   - gives a basic introduction to the Matlab CODE V toolkit
%   cvrgrid    - uses CODE V RAYTRA engine to calculated a grid of rays.
%   cvrgridpol - returns array of polarization ray data at image surface in CodeV
%   cvrmswe    - gets RMS wavefront information from the RMSWE macro
%   cvroot     - changes current working directory in MATLAB to CODE V toolkit directory
%   cvrpol     - returns polarization ray trace data at image surface in CodeV
%   cvrsi      - gets CODE V ray trace data from RSI command
%   cvsave     - saves current lens file under pathfilename
%   cvsd       - returns surface data in a cell array
%   cvsen      - gets CODE V rigid body motion wavefront sensitivity data
%   cvshift    - "shifts" the decenters relative to the current value in CODE V 
%   cvsl       - returns surface number and label in a string
%   cvspot     - returns x,y spot diagram data
%   cvstop     - stops a running CodeV command
%   cvtitle    - queries titles in CODE V
%   cvunits    - gets units of current lens in CodeV
%   cvw        - queries and sets wavelengths and weights in CODE V
%   cvwav      - gets wavefront analysis data from CodeV
%   cvz        - queries and sets active zoom positions in CODE V
