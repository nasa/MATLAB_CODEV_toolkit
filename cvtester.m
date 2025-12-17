% cvtester.m is a script to test the MATLAB_CODEV_toolkit

help 'MATLAB_CODEV_toolkit'

[CODEVpath,defaultspath,CODEVserver] = cvsetup()
cvcv       % starts GUI session of CODE V from Matlab

CVtoolkitPath = fileparts(which('cvon.m'))
cd(CVtoolkitPath)
dir

%start background session of CODE V
cvoff      % kills any running CODE V COM links
cvon       % starts the COM link between Matlab and CODE V

% Load sample lens
cvcmd('res cv_lens:dbgauss')
cvview     %- draws the current lens in Matlab

% Buffer exchange commands
a = flipud(eye(5)); %this is simple check to show up-down difference between matlab and codev
buffer = 1;
cvbufput(a,buffer)
cvbuflis(buffer)
b = cvbuf(1)
b(1,:)
cvcmd('buf put b2 "test"')
cvbuftyp(2)                

cvbufput(pi,buffer)
cvbuf(buffer)
cvdb('num s')
cveva('(num s)')


% Genral help routines and utilities
cvhelp                  
cvlicense               
cvpath                  
cvreadme                
cvroot                  
cvout                   
cvdir                   
cvsave('test.seq')                  
cvstop                  

% System data 
f = cvf        %- retrieves current field points, sets new field points in CodeV,
w = cvw        %- queries and sets wavelengths and weights in CODE V
z = cvz        %- queries and sets active zoom positions in CODE V
cvtitle                 
cvsur                    
cvgc                    
cvdec                   
cvsl                    
cvnum                   
cvunits                 

% Raytrace analysis routines
cvview                  
cvbestsph             
cvenc                   
cvmap                   
cvopl                   
cvpin                   
cvpma                   
cvpsf                   
cvr                     
cvrac                   
cvrbshift               
cvrgrid                 
cvrmswe                 
cvrpol                  
cvrsi                   
cvsen                   
cvshift                 
cvspot                  
