MATLAB-Code V Toolkit
(GSC-15140-1)
 
Overview: 
This toolkit is a set of MATLAB scripts and functions that enable rapid transfer of optical system and performance data from Code V optical software into the MATLAB environment. Typical applications include: extracting prescription data into MATLAB to confirm consistency of various delivered models; perturbing the models and performing various analyses such as ray tracing or generation of point-spread functions in support of integrated modeling activities; and enabling a MATLAB-driven optical model for integrated system-level modeling of wavefront sensing and control.

Instructions
1) Install CODE V
2) Install MATLAB
3) Start MATLAB, and change your directory to the folder with this local repository
4) Run "cvsetup(1)" on the MATLAB command line, and follow instructions to pick your CODE V version and CVUSER folder location.  Note, your "defaults.seq" file is assumed to be in the CVUSER folder, if you have created one.
5) You can now start a background session of CODE V by typing "cvon" in the MATLAB command window to start the process, then "cvin" to load a lens.
6) In general it is good practice to close the session of CODE V by typing "cvoff".

Sample session of commands in Matlab

cvon    %starts CODE V background session
cvin    %opens a gui to locate the lens you want to open
cvnum   %returns the number of surfaces, fields, etc.
cvr     %executes a raytrace
