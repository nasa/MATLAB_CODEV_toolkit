function [data,y,z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP] = cvr(f,r,s,w,z,g,out_option,op_s,mrad)
%CVR  gets CODE V ray trace data from database
%
% FOR REFERENCE RAY DATA ONLY (No relative pupil or field input allowed)
%   function [data,y,z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP] = cvr(f,r,s,w,z,g,out_option,op_s,mrad)
%   
%   INPUTS:  f = field, default = 1, 0 = all fields
%               NOTE:  multiple fields should be input as a single column
%               of fields, since a row is interpreted as the ALT INPUT
%               discussed below.
%            r = ray, 1,2,3,4, or 5:  chief,top,bottom,left,right, default = 1
%                ALT ref-rays:  6,7,8,9:  top right, bottom right, bottom left, top left
%                             r6 = sqrt(2) * [1 1], r7 = sqrt(2) * [1 -1], etc. 
%              = 0, sets r=1:5 
%              = -1, sets r=1:6 
%              = -2, sets r=1:9 
%              NOTE: using r>5 sets the "reference rays" in CODE V 
%            s = surface, default = sI (image), can be multiple surfaces (vector input)
%                use <0 to apply all surfaces
%                can input surface label if single surface desired, e.g. "PRIMARY"
%            w = wave, default = 1, 0 = all waves, -1 = reference wavelength
%            z = zoom, default = 1, 0 = all zooms
%            g = global reference for data, default = -1 (LOCAL)
%            out_option = 0, single structure output, with fields for each raydata
%                            component with size[ f r s w z ]
%                       = 1, default nesting: [raydata f r s w z]
%                       = 2, multi-output, squeezed default nesting
%                       = 3, re-ordered with f w z at end:  [raydata r s f w z]
%                       = 4, xy data only  [xy r s f w z]
%                       = 5, (f,w,z) structure output with fields for each raydata
%                            component with size [s r]
%                       = 6, single structure output, ss rr ff ww zz
%                       = 7, re-ordered with f w z s at end, useful for single surface: [raydata r f w z]
%            op_s = [s1 s2], surfaces to calculate optical path between
%                 default = [1 nums]
%            mrad = 0, output in lens units (default)
%                 = 1, output in meters / radians
% 
%   OUTPUT: 
%           raydata = X,Y,Z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP
%                   nested in an array:
%                       data = [raydata (r)(s)(f)(w)(z)]
%                 NOTE: previous version data order was:
%                       data = shiftdim(raydata(f)(r)(s)(w)(z))
%                 NOTE: OP is given from s1..s except when s=0, where OP is
%                       given from s0..s1.  (This allows systems with the
%                       object at infinity to give back a value.
%           y,z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP:  
%                   >1 output parses data into the requested fields
%                       data is reduced to x
%                       e.g.   [x,y] = cvr(f,r,s,w,z,g);
%   USAGE NOTE:
%       RSI output flags:   A = outside default aperture
%                           E = edge thickness error
%                           O = obstructed
% CODE V Help:
% If the requested ray fails to trace, a diagnostic is given. There are
% three possibilities for ray failure: 
% 1.  The ray encountered a total
% internal reflection (or refracted when TIRO mode was specified).
% 2.  The ray missed a surface.
% 3.  The iteration to a special surface failed to converge.
% 
% For non-sequential surfaces there are also other ray failure modes. In
% any of these cases, the ray is traced up to the failure point. In the
% case of an RSI, if the chief ray for the requested relative field fails
% to trace, then a message to that effect is given, and there is no ray
% trace output. If a ray fails to trace for any reason, a Macro-PLUS
% database item RER (ray error) is set to the surface number of the ray
% failure; the sign tells whether it was the requested ray ( + ) or the
% chief ray ( - ) that failed. If the ray traces correctly, RER is set to
% zero. Three other diagnostics are supplied with the ray trace output: 
% 
% •  The first one relates to negative edge thicknesses. If the ray length
% from the last surface is negative at any surface in the lens, indicating
% virtual ray tracing (tracing backwards along the ray), then an E (for
% edge error) is listed next to the surface number in the ray trace output.
% The usual reason for this is negative edge thickness on a lens, but it
% also can occur in decentered systems, especially on dummy surfaces used
% to reorient the coordinate system (in which case the E flag can be
% ignored).
% 
% •  The other two diagnostics relate to the ray encountering an aperture
% or an obscuration. The output of these diagnostics depends on the setting
% of the Apertures Used field in the System Data window (Lens menu), or CA
% command (CA NO, CA YES, or CA APE). In CODE V a ray traced with either
% RSI or SIN is traced to the image plane without regard to apertures or
% obscurations. In the ray trace output, however, if a ray hits a surface
% outside the apertures (user-defined or default) and the CA status is YES
% (the default for CA), an A (for aperture) is listed next to the surface
% number. If the ray encounters an obscuration, an O (for obscuration) is
% listed next to the surface number. Note that the E diagnostic will
% override the A or O diagnostics in the output. If the CA status is APE
% (for user-defined apertures only), then the check is only made on
% user-defined apertures, i.e., default apertures are ignored. If the CA
% status is NO, then no check is made at all.
% 
% If the ray encounters an aperture or obscuration and sets the A or O
% diagnostic, a Macro-PLUS database item called BLS (blocking surface) is
% set to the surface number. If the ray passes outside the aperture on more
% than one surface, BLS is set to the number of the first such surface. In
% the case of the ray hitting a surface inside an obscuration, BLS is set
% to the negative of the surface number. If the ray encounters neither
% (i.e., passes inside all apertures and outside all obscurations), BLS is
% set to zero. Thus, macros can interrogate RER and BLS to determine the
% status of traced rays.

[nums,numf,numw,numz] = cvnum();
if nargin<1, f = 1; end
if nargin<2, r = 1; end
if nargin<3, s = nums; end
if nargin<4, w = 1; end
if nargin<5, z = 1; end
if nargin<6, g = -1; end
if nargin<7, out_option = 1; end
if nargin<8, op_s = [1 nums]; end
if nargin<9, mrad = false; end
   
% convert empty inputs to defaults
if isempty(f), f = 1; end
if isempty(r), r = 1; end
if isempty(s), s = nums; end
if isempty(w), w = 1; end
if isempty(z), z = 1; end
if isempty(g), g = -1; end
if isempty(out_option), out_option = 1; end
if isempty(op_s) || length(op_s)<2 , op_s = [1 nums]; end

%lookup s if surface name is given
s = convertStringsToChars(s);
if ischar(s), s = cvsl(s,-1); end

%convert 0 or <0 inputs to ALL
if f<1, f = 1:numf; end 
if r==0, r = 1:5;
elseif r==-1, r = 1:6;
elseif r==-2, r = 1:9; 
end
if s<0, s = 0:nums; end
if w==0, w = 1:numw; elseif w==-1, w = cveva('(ref)'); end
if z==0, z = 1:numz; end

CVinput = zeros( max(length(s),max(length(f),max(length(w),length(z)))) , 5 ); %format for cvr.seq
for i=1:length(f), CVinput(i,1) = f(i); end
for i=1:length(r), CVinput(i,2) = r(i); end
for i=1:length(s), CVinput(i,3) = s(i); end
for i=1:length(w), CVinput(i,4) = w(i); end
for i=1:length(z), CVinput(i,5) = z(i); end

%setup reference rays if r>5
if any(r>5)
    
end

cvbufput(CVinput,8);
rows1 = length(f)*length(r)*length(s)*length(w)*length(z);
cvmacro('cvr.seq', [length(f) length(r) length(s) length(w) length(z) g op_s(:)']);

output = cvbuf(1,1:rows1,1:12);
% data = reshape(output',12,length(f),length(r),length(s),length(w),length(z));
data = reshape(output,length(z),length(w),length(s),length(r),length(f),12); %cvr.seq loop structure reverse order
data = permute(data,[6 5 4 3 2 1]); %data is now in order size: raydata f r s w z

% ii=1;
% data = zeros(12,length(f),length(r),length(s),length(w),length(z));
% for ff=1:length(f)
%     for rr=1:length(r) 
%         for ss=1:length(s) 
%             for ww=1:length(w) 
%                 for zz=1:length(z)
%                     data(:,ff,rr,ss,ww,zz) = output(ii,:); %Max is a 6D array containing all data
%                     ii=ii+1;
%                 end
%             end
%         end
%     end
% end

[lensunits2mm,lensunits_str] = cvunits;
if mrad  %convert to meters/radians
    data([1:3 12],:,:,:,:,:) = 1e-3*lensunits2mm*data([1:3 12],:,:,:,:,:); % x y z OP
    data(7:8,:,:,:,:,:) = pi/180*data(7:8,:,:,:,:,:); % AOI AOR
end

if out_option == 0 %single structure output, ff rr ss ww zz
   rdata.datatype   = 'cvr.m'; 
   rdata.field_num   = f; 
   rdata.ray_num   = r; 
   rdata.surface_num   = s; 
   rdata.wave_num   = w; 
   rdata.zoom_num   = z; 
   rdata.global   = g; 
   if mrad
       rdata.units = 'M rads';
   else
       rdata.units = [lensunits_str ' deg'];
   end
   rdata.datasize = '[surfaces rays]';
   rdata.x   = shiftdim(data(1,:,:,:,:,:),1);
   rdata.y   = shiftdim(data(2,:,:,:,:,:),1);
   rdata.z   = shiftdim(data(3,:,:,:,:,:),1);
   rdata.L   = shiftdim(data(4,:,:,:,:,:),1);
   rdata.M   = shiftdim(data(5,:,:,:,:,:),1);
   rdata.N   = shiftdim(data(6,:,:,:,:,:),1);
   rdata.AOI = shiftdim(data(7,:,:,:,:,:),1);
   rdata.AOR = shiftdim(data(8,:,:,:,:,:),1);
   rdata.SRL = shiftdim(data(9,:,:,:,:,:),1);
   rdata.SRM = shiftdim(data(10,:,:,:,:,:),1);
   rdata.SRN = shiftdim(data(11,:,:,:,:,:),1);
   rdata.OP  = shiftdim(data(12,:,:,:,:,:),1);
   if mrad
       rdata.units = 'M rads';
   else
       rdata.units = [lensunits_str ' deg'];
   end
   data = rdata; 

elseif out_option == 6 %single structure output, ss rr ff ww zz
    data = permute(data,[1 4 3 2 5 6]); %data is now in order size: raydata s r f w z
    rdata.datatype   = 'cvr.m';
    rdata.surfaces   = s;
    rdata.rays   = r;
    rdata.fields   = f;
    rdata.waves   = w;
    rdata.zooms   = z;
    rdata.global   = g;
    if mrad
        rdata.units = 'M rads';
    else
        rdata.units = [lensunits_str ' deg'];
    end
    rdata.datasize = '[surfaces rays fields waves zooms]';
    rdata.x   = shiftdim(data(1,:,:,:,:,:),1);
    rdata.y   = shiftdim(data(2,:,:,:,:,:),1);
    rdata.z   = shiftdim(data(3,:,:,:,:,:),1);
    rdata.L   = shiftdim(data(4,:,:,:,:,:),1);
    rdata.M   = shiftdim(data(5,:,:,:,:,:),1);
    rdata.N   = shiftdim(data(6,:,:,:,:,:),1);
    rdata.AOI = shiftdim(data(7,:,:,:,:,:),1);
    rdata.AOR = shiftdim(data(8,:,:,:,:,:),1);
    rdata.SRL = shiftdim(data(9,:,:,:,:,:),1);
    rdata.SRM = shiftdim(data(10,:,:,:,:,:),1);
    rdata.SRN = shiftdim(data(11,:,:,:,:,:),1);
    rdata.OP  = shiftdim(data(12,:,:,:,:,:),1);
    data = rdata;
    
elseif out_option == 5 % (f,w,z) array structure output, rr ss
    for i_f = 1:length(f)
        for i_w = 1:length(w)
            for i_z = 1:length(z)
               rdata(i_f,i_w,i_z).datatype   = 'cvr.m'; 
               rdata(i_f,i_w,i_z).field_num  = f(i_f);
               rdata(i_f,i_w,i_z).surface_num = s;
               rdata(i_f,i_w,i_z).ray_num  = r;
               rdata(i_f,i_w,i_z).wave_num = w(i_w); 
               rdata(i_f,i_w,i_z).zoom_num = z(i_z); 
               rdata(i_f,i_w,i_z).g = g; 
               if mrad
                   rdata(i_f,i_w,i_z).units = 'M rads';
               else
                   rdata(i_f,i_w,i_z).units = [lensunits_str ' deg'];
               end
               rdata(i_f,i_w,i_z).size = '[num_surfaces num_rays]';
               rdata(i_f,i_w,i_z).x   = shiftdim(data(1,i_f,:,:,i_w,i_z),2)'; 
               rdata(i_f,i_w,i_z).y   = shiftdim(data(2,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).z   = shiftdim(data(3,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).L   = shiftdim(data(4,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).M   = shiftdim(data(5,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).N   = shiftdim(data(6,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).AOI = shiftdim(data(7,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).AOR = shiftdim(data(8,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).SRL = shiftdim(data(9,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).SRM = shiftdim(data(10,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).SRN = shiftdim(data(11,i_f,:,:,i_w,i_z),2)';
               rdata(i_f,i_w,i_z).OP  = shiftdim(data(12,i_f,:,:,i_w,i_z),2)';
               if mrad
                   rdata(i_f,i_w,i_z).units = 'M rads';
               else
                   rdata(i_f,i_w,i_z).units = [lensunits_str ' deg'];
               end
            end
        end
    end
    data = rdata;
   
elseif nargout>1  %Multiple output requested (parsed data), ff rr ss ww zz
    
    OP  = shiftdim(data(12,:,:,:,:,:));
    SRN = shiftdim(data(11,:,:,:,:,:));
    SRM = shiftdim(data(10,:,:,:,:,:));
    SRL = shiftdim(data( 9,:,:,:,:,:));
    AOR = shiftdim(data( 8,:,:,:,:,:));
    AOI = shiftdim(data( 7,:,:,:,:,:));
    N   = shiftdim(data( 6,:,:,:,:,:));
    M   = shiftdim(data( 5,:,:,:,:,:));
    L   = shiftdim(data( 4,:,:,:,:,:));
    z   = shiftdim(data( 3,:,:,:,:,:));
    y   = shiftdim(data( 2,:,:,:,:,:));
    x   = shiftdim(data( 1,:,:,:,:,:));

    if out_option==2  %squeeze data output
        OP = squeeze(OP);
        SRN = squeeze(SRN);
        SRM = squeeze(SRM);
        SRL = squeeze(SRL);
        AOR = squeeze(AOR);
        AOI = squeeze(AOI);
        N = squeeze(N);
        M = squeeze(M);
        L = squeeze(L);
        z = squeeze(z);
        y = squeeze(y);
        x = squeeze(x);
    elseif out_option==3 || out_option==4 || out_option==7
        reorder = [2 3 1 4 5];
        OP  = permute(OP ,reorder); % Reorder output to put f,w,z at the end of the data set
        SRN = permute(SRN,reorder); % Reorder output to put f,w,z at the end of the data set
        SRM = permute(SRM,reorder); % Reorder output to put f,w,z at the end of the data set
        SRL = permute(SRL,reorder); % Reorder output to put f,w,z at the end of the data set
        AOR = permute(AOR,reorder); % Reorder output to put f,w,z at the end of the data set
        AOI = permute(AOI,reorder); % Reorder output to put f,w,z at the end of the data set
        N   = permute(N  ,reorder); % Reorder output to put f,w,z at the end of the data set
        M   = permute(M  ,reorder); % Reorder output to put f,w,z at the end of the data set
        L   = permute(L  ,reorder); % Reorder output to put f,w,z at the end of the data set
        z   = permute(z  ,reorder); % Reorder output to put f,w,z at the end of the data set
        y   = permute(y  ,reorder); % Reorder output to put f,w,z at the end of the data set
        x   = permute(x  ,reorder); % Reorder output to put f,w,z at the end of the data set
    end
    
    data = x;
    
else
    if out_option==2 %squeeze data
        data = squeeze(data);
    elseif out_option==3 %put f w z at end, original order = [data f r s w z]
        data = permute(data,[1 3 4 2 5 6]); % Reorder output to put f,w,z at the end of the data set
    elseif out_option==4 %return x y only
        data = permute(data,[1 3 4 2 5 6]); % Reorder output to put f,w,z at the end of the data set
        data = data(1:2,:,:,:,:,:,:);
    elseif out_option==7 %put f w z at end, original order = [data f r s w z]
        data = permute(data,[1 3 2 5 6 4]); % Reorder output to put f,w,z at the end of the data set
    end

end

return
        
% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     

% archive code

if out_option == 0 %structure output
   rdata.x   = data(1,:,:,:,:,:); 
   rdata.y   = data(2,:,:,:,:,:);
   rdata.z   = data(3,:,:,:,:,:);
   rdata.L   = data(4,:,:,:,:,:);
   rdata.M   = data(5,:,:,:,:,:);
   rdata.N   = data(6,:,:,:,:,:);
   rdata.AOI = data(7,:,:,:,:,:);
   rdata.AOR = data(8,:,:,:,:,:);
   rdata.SRL = data(9,:,:,:,:,:);
   rdata.SRM = data(10,:,:,:,:,:);
   rdata.SRN = data(11,:,:,:,:,:);
   rdata.OP  = data(12,:,:,:,:,:);
   data = rdata; 
   
elseif nargout>1  %Multiple output requested (parsed data)
    
    if out_option==3  
        data = permute(data,[4 2 5 6 3 1]); % Reorder output to put f,w,z at the end of the data set
    else
        data = permute(data,[2 3 4 5 6 1]); %put raydata index at end 
    end
    
    OP = data(:,:,:,:,:,12);
    SRN = data(:,:,:,:,:,11);
    SRM = data(:,:,:,:,:,10);
    SRL = data(:,:,:,:,:,9);
    AOR = data(:,:,:,:,:,8);
    AOI = data(:,:,:,:,:,7);
    N = data(:,:,:,:,:,6);
    M = data(:,:,:,:,:,5);
    L = data(:,:,:,:,:,4);
    z = data(:,:,:,:,:,3);
    y = data(:,:,:,:,:,2);
    data = data(:,:,:,:,:,1);

    if out_option==2  %squeeze data output
        OP = squeeze(OP);
        SRN = squeeze(SRN);
        SRM = squeeze(SRM);
        SRL = squeeze(SRL);
        AOR = squeeze(AOR);
        AOI = squeeze(AOI);
        N = squeeze(N);
        M = squeeze(M);
        L = squeeze(L);
        z = squeeze(z);
        y = squeeze(y);
        data = squeeze(data);
    end
    
else
    if out_option==2 %squeeze data
        data = squeeze(data);
    elseif out_option==3 %put f w z at end
        data = permute(data,[1 3 4 2 5 6]); % Reorder output to put f,w,z at the end of the data set
    elseif out_option==4 %return x y only
        data = data(1:2,:,:,:,:,:,:);
    end

end

return

% test script
cvcmd('res CassRC')
f = 1; r = 1:5; s = -1;,w = 1; z = 1; g = -1; outoption = 6; op_s = []; mrad = 1;
rdata_cvr = cvr(f,r,s,w,z,g,outoption,op_s,mrad)
