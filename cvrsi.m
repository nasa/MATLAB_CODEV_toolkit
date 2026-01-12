function [data,y,z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP,BLS,RER] = cvrsi(f,r,s,w,z,g,out_option,op_s,mrad,s_ref)
%CVRSI gets CODE V RSI ray trace data from database
%
%   function [data,y,z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP] = cvrsi(f,r,s,w,z,g,out_option,op_s,mrad,s_ref)
%   
%   INPUTS:  f = field, 1 per ROW
%              = 1 (default), for field 1
%              < 0, all fields:  (1:numf)' (NOTE THE ROW FORMAT!!) 
%              = [fx fy] = relative field input
%            r = ray, 1 per ROW
%              = 1 (default), for reference ray 1, chief ray
%              < 0, all ref rays:  (1:5)' (NOTE THE ROW FORMAT!!) 
%              = 1,2,3,4, or 5:  chief,top,bottom,left,right, default = 1
%              = 6,7,8, or 9:  Rim rays, NE SE SW NW order
%              = [rx ry],  pupil surface relative coordinate
%              = [rx ry s_ref], reference surface actual coordinate with s_ref input
%                USAGE NOTE: all rays are converted to relative x y pupil
%                coords, so VIGNETTING FACTORS ARE NOT USED.  Use cvr.m to 
%                return ray data with VUY, VLY, VUX, VLX implemented.
%            s = surface, default = sI (image), can be multiple surfaces (vector input)
%                use <0 to apply all surfaces
%            w = wave, default = 1, 0 = all waves
%            z = zoom, default = 1, 0 = all zooms
%            g = global reference for data, default = -1 (LOCAL)
%            out_option = 0, structure output
%                       = 1, default nesting: [raydata f r s w z]
%                       = 2, multi-output, squeezed default nesting
%                       = 3, re-orderd with f w z at end:  [raydata r s f w z]
%                       = 4, xy data only  [xy r s f w z]
%                       = 7, re-ordered with f w z s at end, useful for single surface: [raydata r f w z s]
%            op_s = [s1 s2], surfaces to calculate optical path between
%                 default = [1 nums]
%            mrad = 0, output in lens units
%                 = 1, output in meters / radians
%            s_ref = reference surface for relative ray coords
%                    default = [], not used (i.e. RSI assumes STOP)
%   OUTPUT: 
%           data = X,Y,Z,L,M,N,AOI,AOR,SRL,SRM,SRN,OP,BLS,RER
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
%           BLS returns surface number at which ray is first blocked by an aperture or obscuration
%               <0 means Ray is blocked by obscuration on surface
%               >0 means Ray is blocked by aperture on surface
%               0 = Ray is not blocked
%           RER returns ray error flag
%               0 means Trace successful
%               >0 means Value is surface number where failure occurred 
%               <0 means Ray failed when corresponding chief ray (RSI only)
%               failed to trace through the center of the stop; 
%               |RER| is surface number where failure occurred
%   USAGE NOTE:  Vignetting is ignored with RSI, so use cvr if vignetted ray data is requested.
%
%   See also:  cvr


[nums,numf,numw,numz] = cvnum();
if nargin<1, f = 1; end
if nargin<2, r = 1; end
if nargin<3, s = nums; end
if nargin<4, w = 1; end
if nargin<5, z = 1; end
if nargin<6, g = -1; end
if nargin<7, out_option = 1; end
if nargin<8, op_s = [1 nums]; end
if nargin<9, mrad = 0; end
if nargin<10, s_ref = []; end
   
% convert empty inputs to defaults
if isempty(f), f = 1; end
if isempty(r), r = 1; end
if isempty(s), s = nums; end
if isempty(w), w = 1; end
if isempty(z), z = 1; end
if isempty(g), g = -1; end
if isempty(out_option), out_option = 1; end
if isempty(op_s) || length(op_s)<2 , op_s = [1 nums]; end
if isempty(s_ref), s_ref = -ones(size(r,1),1); end

%convert 0 or <0 inputs to ALL
if isscalar(f) && f<1, f = 1:numf; end
if isscalar(r) && r<1, r = 1:5; end
if s<0, s = 0:nums; end
if w==0, w = 1:numw; end
if z==0, z = 1:numz; end

% fix column format for obvious input errors (i.e. make row into column)
if size(f,2)>2 && size(f,1)==1, f = f'; end
if size(r,2)>3 && size(r,1)==1, r = r'; end

%add second column for f and r "relative" values if not present
novalue = -999;
if size(f,2)<2, f = [f f*0]; f(:,2) = novalue; end
if size(r,2)<2, r = [r r*0]; r(:,2) = novalue; end

%reset number of f r w s z values for loops following
numf = size(f,1);
numr = size(r,1);
nums = length(s);
numw = length(w);
numz = length(z);
nums_ref = length(s_ref);

%create a "ring" of rays if r>9
if size(r,1)<2 && size(r,2)>1 && r(1,2) == novalue
    if r(1,1)>9 
        numr = r(1,1);
        for i=1:numr
            theta = 0:2*pi/numr:2*pi*(numr-1)/numr;
            rad = 0*theta + 1;
            [rx,ry] = pol2cart(theta,rad);
            r = [rx(:) ry(:)];
        end
    end
end

%replace ray numbers 6:9 with relative pupil coords NE SE SW NW
for i=1:numr
    if r(i,2)==novalue
        if r(i,1)==1, r(i,:) = [ 0  0]; end
        if r(i,1)==2, r(i,:) = [ 0  1]; end
        if r(i,1)==3, r(i,:) = [ 0 -1]; end
        if r(i,1)==4, r(i,:) = [ 1  0]; end
        if r(i,1)==5, r(i,:) = [-1  0]; end
        if r(i,1)==6, r(i,:) = [ 1  1]/sqrt(2); end
        if r(i,1)==7, r(i,:) = [ 1 -1]/sqrt(2); end
        if r(i,1)==8, r(i,:) = [-1 -1]/sqrt(2); end
        if r(i,1)==9, r(i,:) = [-1  1]/sqrt(2); end
    end
end

%ensure s_ref is the same size as r input
if size(r,2)==3
    s_ref = r(:,3);
    r = r(:,1:2);
end
while length(s_ref)<size(r,1)+1
    s_ref(length(s_ref)+1) = s_ref(length(s_ref));
end

%RSI input for cvrsi.seq
CVinput = zeros( max(size(r,1),max(length(s),max(size(f,1),max(length(w),length(z))))) , 8 );
for i=1:numf, CVinput(i,1:2) = f(i,1:2); end
for i=1:numr, CVinput(i,3:4) = r(i,1:2); CVinput(i,8) = s_ref(i); end
for i=1:nums, CVinput(i,5) = s(i); end
for i=1:numw, CVinput(i,6) = w(i); end
for i=1:numz, CVinput(i,7) = z(i); end

cvbufput(CVinput,8);

rows1 = numf*numr*nums*numw*numz;
out = cvmacro('cvrsi.seq', [numf numr nums numw numz g op_s(:)']);
if contains(out,'ERROR'), disp(out); return; end

% bufdata = cvbuf(1,1:rows1,1:12);
bufdata = cvbuf(1,1:rows1,1:14); %includes BLS and RER

ii=1;
% data = zeros(12,numf,numr,nums,numw,numz);
data = zeros(14,numf,numr,nums,numw,numz); %includes BLS and RER flags
for ff=1:numf
    for rr=1:numr 
        for ss=1:nums 
            for ww=1:numw 
                for zz=1:numz
                    data(:,ff,rr,ss,ww,zz) = bufdata(ii,:); %Max is a 6D array containing all data
                    ii=ii+1;
                end 
            end 
        end 
    end 
end

[lensunits2mm,lensunits_str] = cvunits;

if mrad  %convert to meters/radians
    data([1:3 12],:,:,:,:,:) = 1e-3*lensunits2mm*data([1:3 12],:,:,:,:,:); % x y z OP
    data(7:8,:,:,:,:,:) = pi/180*data(7:8,:,:,:,:,:); % AOI AOR
end

if out_option == 0 %single structure output, ff rr ss ww zz
   rdata.datatype   = 'cvrsi.m'; 
   rdata.field_def   = f; 
   rdata.ray_def   = r; 
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
   rdata.BLS  = shiftdim(data(13,:,:,:,:,:),1);
   rdata.RER  = shiftdim(data(14,:,:,:,:,:),1);
   if mrad
       rdata.units = 'M rads';
   else
       rdata.units = [lensunits_str ' deg'];
   end
   data = rdata; 

elseif out_option == 6 %single structure output, ss rr ff ww zz
    data = permute(data,[1 4 3 2 5 6]); %data is now in order size: raydata s r f w z
   rdata.datatype   = 'cvrsi.m'; 
   rdata.field_def   = f; 
   rdata.ray_def   = r; 
   rdata.surface_num   = s; 
   rdata.wave_num   = w; 
   rdata.zoom_num   = z; 
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
    rdata.BLS  = shiftdim(data(13,:,:,:,:,:),1);
    rdata.RER  = shiftdim(data(14,:,:,:,:,:),1);
    data = rdata;

%     rdata.x   = data(1,:,:,:,:,:);
%     rdata.y   = data(2,:,:,:,:,:);
%     rdata.z   = data(3,:,:,:,:,:);
%     rdata.L   = data(4,:,:,:,:,:);
%     rdata.M   = data(5,:,:,:,:,:);
%     rdata.N   = data(6,:,:,:,:,:);
%     rdata.AOI = data(7,:,:,:,:,:);
%     rdata.AOR = data(8,:,:,:,:,:);
%     rdata.SRL = data(9,:,:,:,:,:);
%     rdata.SRM = data(10,:,:,:,:,:);
%     rdata.SRN = data(11,:,:,:,:,:);
%     rdata.OP  = data(12,:,:,:,:,:);
%     data = rdata;
   
elseif nargout>1  %Multiple output requested (parsed data), ff rr ss ww zz
    
    RER  = shiftdim(data(14,:,:,:,:,:));
    BLS  = shiftdim(data(13,:,:,:,:,:));
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
        RER = squeeze(RER);
        BLS = squeeze(BLS);
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
    elseif out_option==3  
        reorder = [2 3 1 4 5];
        RER = permute(RER ,reorder); % Reorder output to put f,w,z at the end of the data set
        BLS = permute(BLS ,reorder); % Reorder output to put f,w,z at the end of the data set
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
    elseif out_option==3 %put f w z at end
        data = permute(data,[1 3 4 2 5 6]); % Reorder output to put f,w,z at the end of the data set
    elseif out_option==4 %return x y only
        data = permute(data,[1 3 4 2 5 6]); % Reorder output to put f,w,z at the end of the data set
        data = data(1:2,:,:,:,:,:,:);
    elseif out_option==7 %put f w z at end, original order = [data f r s w z]
        data = permute(data,[1 3 2 5 6 4]); % Reorder output to put f,w,z at the end of the data set
    end

end


% !CODE V SYNTAX HELP
% !RSI     [Sk|Si..j] [Zk] [Wn] x_fract_pupil  y_fract_pupil  x_fract_max_x_field  y_fract_max_y_field
% !      Single trace, relative pupil coordinates, relative field
% !RSI     [Sk|Si..j] [Zk] [Wn] Fm  x_fract_pupil  y_fract_pupil
% !      Single trace, numbered field, relative pupil coordinates
% !RSI     [Sk|Si..j] [Zk] [Wn]  Fm  Rs
% !      Single trace, numbered field, numbered reference ray
% !RSI     [Sk|Si..j] [Zk] [Wn] x_ref_surf y_ref_surf x_fract_max_x_field y_fract_max_y_field num_of_ref_surf
% !      Single trace, actual reference surface coordinates and
% !      relative field; num_of_ref_surf defines temporary
% !      reference surface for trace
% !RSI     [Sk|Si..j] [Zk] [Wn] Fm x_ref_surf y_ref_surf num_of_ref_surf
% !      Single trace, numbered field, actual reference surface
% !      coordinates; num_of_ref_surf defines temporary reference surface for trace
% !RSI     [Sk|Si..j]
% !      Trace ray defined by last RSI command

% CODE V HELP
% What You Need to Know About the Relative Single Ray Trace (RSI Command)
% In the Real Ray Trace dialog box, in the Pupil/Field Ray Trace area, the
% arguments for the relative single ray trace (RSI command) are generally
% the X and Y relative pupil coordinates, and the X and Y relative field
% angles, but short forms can be used for pupil and field. If the desired
% field is one of the 25 fields (F2, F2, F3, F4, F5, etc.) in the
% specification data, the field can be specified in Field number (F
% qualifier). A Reference ray number value (R qualifier) can be used for
% the pupil coordinates to select one of the five reference rays (R1, R2,
% R3, R4, R5) but only if an F qualifier has also been used. Other
% qualifiers for the RSI command are Wavelength (W), Zoom Position (Z), and
% Start Surface and End Surface (S). Rays are traced in the reference
% wavelength unless specified with the W qualifier. The chief ray in the
% reference wavelength for that field is traced first; in addition, if the
% requested wavelength is a wavelength other than the reference wavelength,
% a chief ray is also traced at the requested wavelength; for RSI, these
% chief rays must be traceable. The requested ray is traced for all zoom
% positions; a Z qualifier will limit the output to only one zoom position.
% The default is to list ray trace output for all surfaces from object to
% image, unless limited with an S qualifier. On the command line, an RSI
% command with no qualifiers or data will retrace the last ray traced. If
% only the S qualifier is given, the last ray is traced again, with data
% output for the new surface range. This is helpful when you are making
% changes to the lens data and testing with a ray. If the R and F
% qualifiers are not used, pupil and/or field coordinates are given in
% relative values. The X and Y pupil coordinates are relative to the
% circular pupil for the requested relative field angle. The pupil used is
% the paraxial entrance pupil diameter located at the real entrance pupil
% location for that field. A value of 1 for relative X or relative Y means
% trace the ray at the edge of the pupil in that direction. Vignetting
% values are ignored by RSI. Since the vignetting values change for
% different fields, RSI can not know what the vignetting value would be for
% an arbitrary field, hence RSI does not use any vignetting values. You
% must manually take them into account. This is true even when using the F
% qualifier, even though the specified fields do have defined vignetting
% values. Note also that RSI ignores any CRA specification unless the field
% is specified with an F qualifier, in which case the CRA specification is
% used.
%  
% YAN 0 10 VUY 0 .3 RSI 0 1 0 1 ! Traces the ray at the upper rim of the
% pupil at
%   ! the maximum field (vignetting value is ignored)
% RSI 0 .7 0 1 ! Traces a ray at a relative pupil position
%   ! of 0.7 in Y (the same as the vignetted ! marginal ray) at the full
%   field
% RSI 0 1-(VUY F2) 0 1 ! Traces the upper Y vignetted ray (same as R2)
%   ! for full field
% RSI F2 0 1 ! Traces the ray at the upper rim of the pupil
%   ! (vignetting is ignored) for field 2
% RSI F2 0 1-(VUY F2) ! Traces the upper Y vignetted ray for field 2
%   ! (same as R2)
% 
% If the R qualifier is used, then the reference rays are traced, and these
% reference rays do take into account any vignetting.
%  
% YAN 0 10 VUY 0 .3 RSI F2 R2 ! Traces the upper vignetted marginal ray for
% field 2
%   ! (reference ray R2)
% 
% The field angle inputs to the RSI command are relative to the maximum X
% and Y values for any of the 25 field angles specified in the LDM
% (regardless of whether the largest value is the 25th field or not); the X
% and Y maximums are selected separately and may thus come from different
% field points. The sign of the maximum value is retained as the definition
% of the “largest” field. Thus, if the maximum field angle specified in the
% LDM is negative, then a relative field of +1 in RSI corresponds to the
% negative angle. If there are both a positive and a negative field angle
% with the same absolute maximum value, the positive value is retained for
% the maximum. X and Y relative angles are calculated separately for the X
% and Y field specifications. If there are no X field specifications, then
% any relative X field in a RSI command will trace a 0 field, since the
% maximum value is 0; if you want to trace a ray from an X field, you must
% first have a non-zero X field specification entered in the LDM.
%  
%  
% YAN 0 10 20 RSI 0 0 0 1 ! Traces a ray at +20 degrees
%     
% YAN 0 -10 -20 RSI 0 0 0 1 ! Traces a ray at -20 degrees RSI 0 0 0 -1 !
% Traces a ray at +20 degrees
%     
% YAN 0 -10 10 RSI 0 0 0 1 ! Traces a ray at +10 degrees
%     
% YAN 0 10 XAN 0 0 RSI 0 0 1 0 ! Traces an on-axis ray (no X field
% specified)
%     
% YAN 0 10 XAN 0 5 RSI 0 0 1 0 ! Traces a ray at +5 degrees in X RSI 0 0 1
% 1 ! Traces a ray at 5 degrees in X and 10 degrees in Y
% 
% Since object heights are stored in CODE V, relative angles do not
% multiply the angle value, but multiply the tangent of the maximum angle.
%  
%  
% YAN 0 10 RSI 0 0 0 1 ! Traces a ray at 10 degrees
%   
% RSI 0 0 0 .5 ! Traces a ray at 5.038369 degrees (the angle
%   ! whose tangent is half the tangent of 10 degrees)
% 
% What is used by RSI is the maximum value of the field specifications. It
% does not matter what the other values are or how many there are.
%  
%  
% YAN 10 RSI 0 0 0 .5 ! Trace a ray at 5.038369 degrees YAN 1 2 10 3 4 RSI
% 0 0 0 .5 ! Trace a ray at 5.038369 degrees
% 
% If the F qualifier is used, then both the X and Y field specifications
% are used for that field number. For example, set up three fields on-axis
% and along the diagonal of a square format:
%  
% YAN 0 5 10 XAN 0 5 10
% 
% There is no F qualifier for these fields which will represent a ray at
% 10× in Y only. To trace such a ray, relative fields must be used, as in:
%  
% RSI 0 0 0 1 ! Trace a ray at 10 degrees in Y and 0 degrees in X
% 
% Relative Real Ray Trace – Targeted to a Spot on a Surface The relative
% real ray trace (RSI command) also allows targeting of a ray to a specific
% spot on a specific surface. In the GUI, select the Aim at Surface button
% to enable this option. On the command line, this is accomplished by
% adding a fifth entry to the RSI data arguments. This fifth entry is the
% surface number of the targeted surface (Reference surface field in the
% GUI). If this is used, then the first two entries are not relative pupil
% coordinates, but are actual X and Y positions on the targeted surface
% (Actual coordinates on reference surface in the GUI). CODE V will then
% iterate to find the ray from the specified field which hits the targeted
% surface at the targeted spot. The field may be specified with a Field
% number (F) qualifier. If the F qualifier is used, then the surface number
% is the third entry.
%  
% RSI .3 .5 0 1 7 ! Trace a ray at full field in Y which hits
%   ! surface 7 at an X value of .3 and a Y value of .5
% RSI F1 1 2 5 ! Trace a ray at the first specified field to
%   ! hit surface 5 at a point of 1 in X and 2 in Y
% 
% If the F qualifier is used, then any CRA specifications are used,
% otherwise they are ignored. The CRA will only affect the chief ray
% traced, and will have no effect on the targeted ray. Note that since RSI
% requires the chief ray to trace, the CRA and F qualifier can be used to
% trace a ray that would normally not be traced due to chief ray failure.
% In a non-sequential system, if a ray is requested to hit a surface at a
% specific point, and the ray can hit the surface two or more times, the
% target is applied to the first hit at the surface. (Note that in
% Automatic Design, ray data for the last hit on the surface is the data
% retained and used for any constraints.)


% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     